#!/bin/bash
#
# master_job.sh — Run label-transfer benchmark via targets pipeline (HPC/local)
#
# Orchestrates the label_transfer_task/_targets.R targets pipeline with proper
# renv loading, logging and error handling.
#
# Usage:
#   bash label_transfer_task/scripts/master_job.sh          # Run locally
#   sbatch label_transfer_task/scripts/submit_hpc.sh        # Submit to HPC via SLURM
#
# Requirements:
#   - Run from project root: cd /path/to/ISC_benchmark_reproducibility
#   - R with packages: targets, yaml, dplyr, Seurat, BiocParallel
#   - Processed data in data/processed/ (run data_processing pipeline first)
#

set -euo pipefail

# Configuration
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

log_msg() {
  local timestamp
  timestamp=$(date '+%Y-%m-%d %H:%M:%S')
  echo "[${timestamp}] $1"
}

log_msg "============================================"
log_msg "Starting label-transfer benchmark"
log_msg "============================================"
log_msg "Project root: ${PROJECT_ROOT}"
log_msg "Pipeline directory: ${LT_DIR}"

# Load modules if available (HPC environments)
# Kept in sync with ISC_benchmark/scripts/submit_hpc.sh: R 4.5.2 is the only renv library
# tree with packages like anndataR that the shared python-pipeline helpers need.
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

# Check prerequisites
if [[ ! -f "${LT_DIR}/config/label_transfer_parameters.yaml" ]]; then
  log_msg "ERROR: label_transfer_task/config/label_transfer_parameters.yaml not found."
  exit 1
fi

if [[ ! -d "${PROJECT_ROOT}/data/processed" ]]; then
  log_msg "ERROR: data/processed/ not found. Run data_processing pipeline first."
  exit 1
fi

# Check required R packages under this project's renv
if ! Rscript - <<'RS' 2>&1; then
project_root <- Sys.getenv("PROJECT_ROOT")
stopifnot(nzchar(project_root))

# Hardcoded to match ISC_benchmark's setup: this cluster's renv library only has a
# linux-rocky-9.8/R-4.5 tree (the R.version-derived path silently no-ops otherwise).
renv_lib <- file.path(project_root, "renv", "library", "linux-rocky-9.8", "R-4.5", "x86_64-pc-linux-gnu")
if (dir.exists(renv_lib)) {
  .libPaths(unique(c(renv_lib, .libPaths())))
}

activate <- file.path(project_root, "renv", "activate.R")
if (file.exists(activate)) source(activate)

stopifnot(requireNamespace("renv", quietly = TRUE))
renv::load(project = project_root)

cat("[INFO] .libPaths():\n")
writeLines(.libPaths())

stopifnot(
  requireNamespace("targets",      quietly = TRUE),
  requireNamespace("yaml",         quietly = TRUE),
  requireNamespace("dplyr",        quietly = TRUE),
  requireNamespace("BiocParallel", quietly = TRUE)
)
RS
  log_msg "ERROR: Required R packages not loadable (after renv::load())."
  log_msg "If this is a new machine, run: R -e 'renv::restore()' from the project root."
  exit 1
fi

# Run targets pipeline
log_msg ""
log_msg "Running targets pipeline (label_transfer_task/_targets.R)..."

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

cat("[INFO] .libPaths():\n")
writeLines(.libPaths())

library(targets)
targets::tar_make(callr_function = NULL)
RS
  log_msg ""
  log_msg "✓ Label-transfer benchmark completed successfully"
  log_msg ""
  log_msg "Results:"
  log_msg "  label_transfer_task/results/aggregated/label_transfer_metrics_aggregated.csv"
else
  log_msg ""
  log_msg "✗ Label-transfer benchmark failed. See logs."
  exit 1
fi
