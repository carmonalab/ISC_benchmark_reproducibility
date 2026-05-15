#!/bin/bash
#
# master_job_between_datasets.sh — Run between-dataset label-transfer benchmark
#
set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
LT_DIR="${PROJECT_ROOT}/label_transfer_task"

export PROJECT_ROOT
export LT_DIR

log_msg() {
  local timestamp
  timestamp=$(date '+%Y-%m-%d %H:%M:%S')
  echo "[${timestamp}] $1"
}

log_msg "============================================"
log_msg "Starting between-dataset label-transfer benchmark"
log_msg "============================================"
log_msg "Project root: ${PROJECT_ROOT}"
log_msg "Pipeline directory: ${LT_DIR}"

if command -v module >/dev/null 2>&1; then
  log_msg "Loading environment modules"
  module purge || true
  module load GCC/12.3.0 || true
  module load R/4.3.2 || true
  module load GLPK/5.0 || true
  module load cairo/1.17.8 || true
  module load freetype/2.13.0 || true
  module load libwebp/1.3.1 || true
fi

if [[ ! -f "${LT_DIR}/config/label_transfer_between_parameters.yaml" ]]; then
  log_msg "ERROR: label_transfer_task/config/label_transfer_between_parameters.yaml not found."
  exit 1
fi

if [[ ! -d "${PROJECT_ROOT}/data/processed" ]]; then
  log_msg "ERROR: data/processed/ not found. Run data_processing pipeline first."
  exit 1
fi

log_msg ""
log_msg "Running targets pipeline (_targets_between_datasets.R)..."

cd "${LT_DIR}"
if Rscript - <<'RS' 2>&1; then
options(repos = c(CRAN = "https://packagemanager.posit.co/cran/2024-01-15"))
project_root <- Sys.getenv("PROJECT_ROOT")
stopifnot(nzchar(project_root))

r_mm <- paste0(R.version$major, ".", sub("\\..*$", "", R.version$minor))
renv_lib <- file.path(project_root, "renv", "library", paste0("R-", r_mm), R.version$platform)
if (dir.exists(renv_lib)) {
  .libPaths(unique(c(renv_lib, .libPaths())))
}

activate <- file.path(project_root, "renv", "activate.R")
if (file.exists(activate)) source(activate)
stopifnot(requireNamespace("renv", quietly = TRUE))
renv::load(project = project_root)

library(targets)
targets::tar_make(
  script = "_targets_between_datasets.R",
  store = "_targets_between",
  callr_function = NULL
)
RS
  log_msg ""
  log_msg "✓ Between-dataset label-transfer benchmark completed successfully"
  log_msg ""
  log_msg "Results:"
  log_msg "  label_transfer_task/results_between/aggregated/label_transfer_between_metrics_aggregated.csv"
else
  log_msg ""
  log_msg "✗ Between-dataset label-transfer benchmark failed. See logs."
  exit 1
fi
