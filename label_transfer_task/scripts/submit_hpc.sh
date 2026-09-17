#!/bin/bash
#SBATCH --job-name=ISC_label_transfer
#SBATCH --partition=shared-cpu # Adjust to your HPC partition
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --output=label_transfer_task/logs/label_transfer_%j.log
#SBATCH --error=label_transfer_task/logs/label_transfer_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=josep.garnicacaparros@unige.ch

# ============================================================================
# SLURM Job Submission: Label-Transfer Benchmark via Targets Pipeline
#
# This script submits the label-transfer benchmark job to HPC via SLURM.
# It runs master_job.sh, which orchestrates the targets-based pipeline.
#
# Before running:
#   1. Ensure processed data exists in data/processed/ (run data_processing first)
#   2. Optionally restrict datasets/replicates in:
#      label_transfer_task/config/label_transfer_parameters.yaml
#
# Submit the job (from project root):
#   sbatch label_transfer_task/scripts/submit_hpc.sh
#
# Monitor:
#   squeue -u $USER
#   tail -f label_transfer_task/logs/label_transfer_<jobid>.log
#
# ============================================================================

set -euo pipefail

# Load required modules (kept in sync with ISC_benchmark/scripts/submit_hpc.sh and
# label_transfer_task/scripts/master_job.sh: R 4.5.2 is the only renv library tree
# with packages the shared python-pipeline helpers need, e.g. anndataR).
module purge
module load GCCcore/10.3.0
module load Python/3.9.5-bare
module load GCC/14.3.0
module load R/4.5.2
module load GLPK/5.0 || true
module load cairo/1.17.8 || true
module load freetype/2.13.0 || true
module load libwebp/1.3.1 || true

# Keep the SCCAF venv's Python runtime compatible after newer GCC/R modules are loaded.
export LD_LIBRARY_PATH="/opt/ebsofts/Python/3.11.5-GCCcore-13.2.0/lib:/opt/ebsofts/libffi/3.3-GCCcore-10.3.0/lib64:${LD_LIBRARY_PATH:-}"

# Activate renv (best-effort; master_job.sh will also activate)
PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
export RENV_PROJECT="${PROJECT_ROOT}"
export RENV_PROJECT_EXPLICIT="${PROJECT_ROOT}"
export RENV_CONFIG_AUTOLOADER_ENABLED="FALSE"
if [[ -f "${PROJECT_ROOT}/renv/activate.R" ]]; then
  Rscript -e "source('${PROJECT_ROOT}/renv/activate.R'); renv::load(project='${PROJECT_ROOT}')" >/dev/null 2>&1 || true
fi

# Create logs directory
mkdir -p label_transfer_task/logs

# Run the master job script
bash label_transfer_task/scripts/master_job.sh
