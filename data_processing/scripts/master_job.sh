#!/bin/bash
#
# master_job.sh — Process all datasets via SLURM/HPC
#
# Usage:
#   bash data_processing/scripts/master_job.sh     # Run locally (serial)
#   sbatch data_processing/scripts/submit_hpc.sh   # Submit to HPC
#
# Requirements:
#   - Run from project root: cd /path/to/ISC_benchmark_reproducibility
#   - R with packages: yaml, dplyr, stringr, Seurat, Matrix, BiocParallel
#

set -euo pipefail

# Configuration
PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
DATA_PROCESSING_DIR="${PROJECT_ROOT}/data_processing"

# Logging
LOG_FILE="${DATA_PROCESSING_DIR}/data_processing.log"
: > "${LOG_FILE}"  # Clear log

log_msg() {
  local msg="$1"
  local timestamp=$(date '+%Y-%m-%d %H:%M:%S')
  echo "[${timestamp}] ${msg}" | tee -a "${LOG_FILE}"
}

log_msg "============================================"
log_msg "Starting dataset processing"
log_msg "============================================"
log_msg "Project root: ${PROJECT_ROOT}"
log_msg "Processing directory: ${DATA_PROCESSING_DIR}"

# Check prerequisites
if [[ ! -f "${DATA_PROCESSING_DIR}/config/core_datasets.yaml" ]]; then
  log_msg "ERROR: data_processing/config/core_datasets.yaml not found."
  exit 1
fi

if [[ ! -f "${DATA_PROCESSING_DIR}/config/processing_parameters.yaml" ]]; then
  log_msg "ERROR: data_processing/config/processing_parameters.yaml not found."
  exit 1
fi

if [[ ! -f "${DATA_PROCESSING_DIR}/process_datasets.R" ]]; then
  log_msg "ERROR: data_processing/process_datasets.R not found."
  exit 1
fi

# Check if raw data exists
RAW_DIR="${PROJECT_ROOT}/data/raw"
if [[ ! -d "${RAW_DIR}" ]] || [[ -z "$(find "${RAW_DIR}" -maxdepth 1 -name '*.rds' 2>/dev/null | head -1)" ]]; then
  log_msg "WARNING: No .rds files found in data/raw/"
  log_msg "Download datasets from Zenodo first (see data_processing/README.md)"
fi

# Run R processing script
log_msg ""
log_msg "Running R processing script..."
log_msg "Command: Rscript ${DATA_PROCESSING_DIR}/process_datasets.R"

cd "${PROJECT_ROOT}"
if Rscript "${DATA_PROCESSING_DIR}/process_datasets.R" 2>&1 | tee -a "${LOG_FILE}"; then
  log_msg ""
  log_msg "✓ Dataset processing completed successfully"
else
  log_msg ""
  log_msg "✗ Dataset processing failed. See log: ${LOG_FILE}"
  exit 1
fi

