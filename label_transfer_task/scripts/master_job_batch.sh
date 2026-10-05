#!/bin/bash
#
# master_job_batch.sh — Run batch label-transfer benchmark
#
set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
LT_DIR="${PROJECT_ROOT}/label_transfer_task"

export PROJECT_ROOT
export LT_DIR
# Pin renv's project detection to PROJECT_ROOT; otherwise `cd "${LT_DIR}"` below makes
# renv/activate.R treat label_transfer_task itself as an uninitialized project and
# bootstrap a brand-new renv/ scaffold there.
export RENV_PROJECT="${PROJECT_ROOT}"
export RENV_PROJECT_EXPLICIT="${PROJECT_ROOT}"
export RENV_CONFIG_AUTOLOADER_ENABLED="FALSE"

# Branches run in parallel via future::multicore (fork-based); forking a process
# with an already-threaded BLAS causes "stack imbalance" warnings and CPU
# oversubscription, so keep BLAS/OMP single-threaded per worker.
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export FLEXIBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

log_msg() {
  local timestamp
  timestamp=$(date '+%Y-%m-%d %H:%M:%S')
  echo "[${timestamp}] $1"
}

log_msg "============================================"
log_msg "Starting batch label-transfer benchmark"
log_msg "============================================"
log_msg "Project root: ${PROJECT_ROOT}"
log_msg "Pipeline directory: ${LT_DIR}"

# Kept in sync with label_transfer_task/scripts/master_job.sh: R 4.5.2 is the only renv
# library tree matching the root renv.lock (e.g. contains anndataR).
if command -v module >/dev/null 2>&1; then
  log_msg "Loading environment modules"
  module purge || true
  module load GCCcore/10.3.0 || true
  module load Python/3.9.5-bare || true
  module load GCC/14.3.0 || true
  module load R/4.5.2 || true
  module load GLPK/5.0 || true
  module load cairo/1.17.8 || true
  module load freetype/2.13.0 || true
  module load libwebp/1.3.1 || true
fi

# Keep the SCCAF venv's Python runtime compatible after newer GCC/R modules are loaded.
export LD_LIBRARY_PATH="/opt/ebsofts/Python/3.9.5-GCCcore-10.3.0/lib:/opt/ebsofts/libffi/3.3-GCCcore-10.3.0/lib64:${LD_LIBRARY_PATH:-}"

if [[ ! -f "${LT_DIR}/config/label_transfer_batch_parameters.yaml" ]]; then
  log_msg "ERROR: label_transfer_task/config/label_transfer_batch_parameters.yaml not found."
  exit 1
fi

if [[ ! -d "${PROJECT_ROOT}/data/processed" ]]; then
  log_msg "ERROR: data/processed/ not found. Run data_processing pipeline first."
  exit 1
fi

log_msg ""
log_msg "Running targets pipeline (_targets_between_batch.R)..."

cd "${LT_DIR}"
if Rscript - <<'RS' 2>&1; then
options(repos = c(CRAN = "https://packagemanager.posit.co/cran/2024-01-15"))
project_root <- Sys.getenv("PROJECT_ROOT")
stopifnot(nzchar(project_root))

renv_lib <- file.path(project_root, "renv", "library", "linux-rocky-9.8", "R-4.5", "x86_64-pc-linux-gnu")
if (dir.exists(renv_lib)) {
  .libPaths(unique(c(renv_lib, .libPaths())))
}

activate <- file.path(project_root, "renv", "activate.R")
if (file.exists(activate)) source(activate)
stopifnot(requireNamespace("renv", quietly = TRUE))
renv::load(project = project_root)

library(targets)
library(future)

n_workers <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "4"))
cat("[INFO] Running branches in parallel with", n_workers, "future workers\n")
future::plan(future::multicore, workers = n_workers)
targets::tar_make_future(
  script = "_targets_between_batch.R",
  store = "_targets_batch",
  workers = n_workers,
  callr_function = NULL
)
RS
  log_msg ""
  log_msg "✓ Between-dataset label-transfer benchmark completed successfully"
  log_msg ""
  log_msg "Results:"
  log_msg "  label_transfer_task/results_batch/aggregated/label_transfer_batch_metrics_aggregated.csv"
else
  log_msg ""
  log_msg "✗ Between-dataset label-transfer benchmark failed. See logs."
  exit 1
fi
