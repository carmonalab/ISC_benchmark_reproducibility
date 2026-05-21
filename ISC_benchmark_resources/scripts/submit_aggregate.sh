#!/bin/bash -l
# Aggregate resource benchmarking results after all jobs complete
# This can be submitted as a dependent job using:
#   sbatch --dependency=afterok:JOB_ID scripts/submit_aggregate.sh

set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
RESOURCES_DIR="${PROJECT_ROOT}/ISC_benchmark_resources"
export PROJECT_ROOT
export RENV_PROJECT="${PROJECT_ROOT}"
export RENV_CONFIG_AUTOLOADER_ENABLED="FALSE"

log_msg() {
  local timestamp
  timestamp=$(date '+%Y-%m-%d %H:%M:%S')
  echo "[${timestamp}] $1"
}

log_msg "============================================"
log_msg "Aggregating ISC resource benchmark results"
log_msg "============================================"
log_msg "Project root: ${PROJECT_ROOT}"

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

cd "${RESOURCES_DIR}"

# Run aggregation script
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

source("R/aggregate_results.R")
RS
  log_msg "✓ Aggregation completed successfully"
else
  log_msg "✗ Aggregation failed. See logs."
  exit 1
fi
