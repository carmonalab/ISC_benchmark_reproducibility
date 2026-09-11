#!/bin/bash -l
#SBATCH --job-name=isc_resources_aggregate
#SBATCH --partition=public-cpu
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --output=ISC_benchmark_resources/logs/aggregate-%j.out
#SBATCH --error=ISC_benchmark_resources/logs/aggregate-%j.err
# Aggregate resource benchmarking results after all jobs complete
# This can be submitted as a dependent job using:
#   sbatch --dependency=afterok:JOB_ID scripts/submit_aggregate.sh

set -euo pipefail

PROJECT_ROOT="${PROJECT_ROOT:-${SLURM_SUBMIT_DIR:-}}"
if [[ ! -d "${PROJECT_ROOT}/ISC_benchmark_resources" ]]; then
  PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
fi
RESOURCES_DIR="${PROJECT_ROOT}/ISC_benchmark_resources"
export PROJECT_ROOT
export RENV_PROJECT="${PROJECT_ROOT}"
export RENV_PROJECT_EXPLICIT="${PROJECT_ROOT}"
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
  module purge 2>/dev/null || true
  module load GCCcore/10.3.0 2>/dev/null || true
  module load Python/3.9.5-bare 2>/dev/null || true
  module load GCC/14.3.0 2>/dev/null || true
  module load R/4.5.2 2>/dev/null || true
  module load GLPK/5.0 2>/dev/null || true
  module load cairo/1.17.8 2>/dev/null || true
  module load freetype/2.13.0 2>/dev/null || true
  module load libwebp/1.3.1 2>/dev/null || true
fi

# Keep Python runtimes for the external tool venvs compatible after newer GCC/R modules are loaded.
export LD_LIBRARY_PATH="/opt/ebsofts/Python/3.11.5-GCCcore-13.2.0/lib:/opt/ebsofts/libffi/3.3-GCCcore-10.3.0/lib64:${LD_LIBRARY_PATH:-}"

cd "${RESOURCES_DIR}"

# Run aggregation script
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

source("R/aggregate_results.R")
RS
  log_msg "✓ Aggregation completed successfully"
else
  log_msg "✗ Aggregation failed. See logs."
  exit 1
fi
