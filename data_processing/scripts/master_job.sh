#!/bin/bash
#
# master_job.sh — Process all datasets via targets pipeline (HPC/local)
#
# This script orchestrates the targets-based data processing pipeline
# (data_processing/_targets.R) with proper logging and error handling.
#
# Usage:
#   bash data_processing/scripts/master_job.sh          # Run locally (serial)
#   sbatch data_processing/scripts/submit_hpc.sh        # Submit to HPC via SLURM
#
# Requirements:
#   - Run from project root: cd /path/to/ISC_benchmark_reproducibility
#   - R with packages: targets, yaml, dplyr, stringr, Seurat, Matrix, BiocParallel
#   - Raw data downloaded to data/raw/ (see data_processing/README.md)
#

set -euo pipefail

# Configuration
PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
DATA_PROCESSING_DIR="${PROJECT_ROOT}/data_processing"

# Export for R subprocesses
export PROJECT_ROOT
export DATA_PROCESSING_DIR

log_msg() {
  local msg="$1"
  local timestamp=$(date '+%Y-%m-%d %H:%M:%S')
  echo "[${timestamp}] ${msg}"
}

log_msg "============================================"
log_msg "Starting dataset processing"
log_msg "============================================"
log_msg "Project root: ${PROJECT_ROOT}"
log_msg "Processing directory: ${DATA_PROCESSING_DIR}"

# Load modules if available (HPC environments)
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

# Check prerequisites
if [[ ! -f "${DATA_PROCESSING_DIR}/config/core_datasets.yaml" ]]; then
  log_msg "ERROR: data_processing/config/core_datasets.yaml not found."
  exit 1
fi

if [[ ! -f "${DATA_PROCESSING_DIR}/config/processing_parameters.yaml" ]]; then
  log_msg "ERROR: data_processing/config/processing_parameters.yaml not found."
  exit 1
fi

# Check if required R packages are available (under this project's renv)
if ! Rscript - <<'RS' 2>&1; then
project_root <- Sys.getenv("PROJECT_ROOT")
stopifnot(nzchar(project_root))

# Force project renv library on .libPaths() (works even before renv::load())
r_mm <- paste0(R.version$major, ".", sub("\\..*$", "", R.version$minor))
renv_lib <- file.path(project_root, "renv", "library", paste0("R-", r_mm), R.version$platform)
if (dir.exists(renv_lib)) {
  .libPaths(unique(c(renv_lib, .libPaths())))
}

# Ensure renv is activated for this project, then load it.
activate <- file.path(project_root, "renv", "activate.R")
if (file.exists(activate)) {
  source(activate)
}

stopifnot(requireNamespace("renv", quietly = TRUE))
renv::load(project = project_root)

cat("[INFO] .libPaths():\n")
writeLines(.libPaths())

stopifnot(
  requireNamespace("targets", quietly = TRUE),
  requireNamespace("yaml", quietly = TRUE),
  requireNamespace("Seurat", quietly = TRUE),
  requireNamespace("BiocParallel", quietly = TRUE)
)
RS
  log_msg "ERROR: Required R packages not loadable in this environment (after renv::load())."
  log_msg "If this is a new machine, run: R -e 'renv::restore()' from the project root."
  exit 1
fi

# Check if raw data exists
RAW_DIR="${PROJECT_ROOT}/data/raw"
if [[ ! -d "${RAW_DIR}" ]] || [[ -z "$(find "${RAW_DIR}" -maxdepth 1 -name '*.rds' 2>/dev/null | head -1)" ]]; then
  log_msg "WARNING: No .rds files found in data/raw/"
  log_msg "Download datasets from Zenodo first (see data_processing/README.md)"
fi

# Run targets pipeline
log_msg ""
log_msg "Running targets pipeline (data_processing/_targets.R)..."
log_msg "Command: cd ${DATA_PROCESSING_DIR} && Rscript -e 'source(\"${PROJECT_ROOT}/renv/activate.R\"); renv::load(project=\"${PROJECT_ROOT}\"); targets::tar_make(callr_function=NULL)'"

cd "${DATA_PROCESSING_DIR}"
if Rscript - <<'RS' 2>&1; then
options(repos = c(CRAN = "https://packagemanager.posit.co/cran/2024-01-15"))
project_root <- Sys.getenv("PROJECT_ROOT")
stopifnot(nzchar(project_root))

# Force project renv library on .libPaths() (works even before renv::load())
r_mm <- paste0(R.version$major, ".", sub("\\..*$", "", R.version$minor))
renv_lib <- file.path(project_root, "renv", "library", paste0("R-", r_mm), R.version$platform)
if (dir.exists(renv_lib)) {
  .libPaths(unique(c(renv_lib, .libPaths())))
}

activate <- file.path(project_root, "renv", "activate.R")
if (file.exists(activate)) {
  source(activate)
}
stopifnot(requireNamespace("renv", quietly = TRUE))
renv::load(project = project_root)

cat("[INFO] .libPaths():\n")
writeLines(.libPaths())
library(targets)
targets::tar_make(callr_function = NULL)

tryCatch({
  cat("\n[SUMMARY]\n")
  print(targets::tar_read(summary))
}, error = function(e) {
  message("[INFO] Could not read summary target: ", conditionMessage(e))
})
RS
  log_msg ""
  log_msg "✓ Dataset processing completed successfully"
  log_msg ""
  log_msg "View results:"
  log_msg "  targets::tar_read(summary)     # Summary of all datasets"
  log_msg "  targets::tar_read(report)      # Detailed processing report"
else
  log_msg ""
  log_msg "✗ Dataset processing failed. See SLURM output/error logs."
  exit 1
fi

