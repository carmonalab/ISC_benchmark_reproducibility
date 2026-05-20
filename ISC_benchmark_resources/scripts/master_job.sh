#!/bin/bash
set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
RESOURCES_DIR="${PROJECT_ROOT}/resources"

export PROJECT_ROOT
export RESOURCES_DIR

log_msg() {
  local timestamp
  timestamp=$(date '+%Y-%m-%d %H:%M:%S')
  echo "[${timestamp}] $1"
}

log_msg "============================================"
log_msg "Starting ISC resource benchmark"
log_msg "============================================"
log_msg "Project root: ${PROJECT_ROOT}"
log_msg "Pipeline directory: ${RESOURCES_DIR}"
if [[ -n "${RESOURCE_DATASET_ID:-}" ]]; then
  log_msg "Dataset filter: ${RESOURCE_DATASET_ID}"
fi

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

if [[ ! -f "${RESOURCES_DIR}/config/resource_parameters.yaml" ]]; then
  log_msg "ERROR: resources/config/resource_parameters.yaml not found."
  exit 1
fi

if [[ ! -d "${PROJECT_ROOT}/data/processed" ]]; then
  log_msg "ERROR: data/processed/ not found. Run data_processing pipeline first."
  exit 1
fi

if ! Rscript - <<'RS' 2>&1; then
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

stopifnot(
  requireNamespace("targets", quietly = TRUE),
  requireNamespace("yaml", quietly = TRUE),
  requireNamespace("dplyr", quietly = TRUE),
  requireNamespace("scTypeEval", quietly = TRUE)
)
RS
  log_msg "ERROR: Required R packages not loadable after renv::load()."
  exit 1
fi

cd "${RESOURCES_DIR}"

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

selector <- Sys.getenv("RESOURCE_DATASET_ID", unset = "")
if (!nzchar(selector)) {
  selector <- Sys.getenv("RESOURCE_DATASET_IDS", unset = "")
}
if (!nzchar(selector)) {
  selector <- "default"
}

safe_selector <- gsub("[^A-Za-z0-9._-]", "_", selector)
store_dir <- file.path("_targets", paste0("store_", safe_selector))
dir.create(store_dir, showWarnings = FALSE, recursive = TRUE)

message("Starting targets::tar_make(callr_function = NULL) at ", Sys.time())
targets::tar_make(callr_function = NULL, store = store_dir)
message("Finished targets::tar_make(callr_function = NULL) at ", Sys.time())
RS
  log_msg "✓ ISC resource benchmark completed successfully"
  log_msg "Outputs: resources/output/<dataset>/<ident>/*.rds"
else
  log_msg "✗ ISC resource benchmark failed. See logs."
  exit 1
fi