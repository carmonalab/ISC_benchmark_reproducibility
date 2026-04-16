#!/bin/bash
#SBATCH --job-name=ISC_label_transfer
#SBATCH --partition=private-carmona-gpu
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=256G
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

# Load required modules (customize for your HPC)
module purge
module load GCC/12.3.0
module load R/4.3.2
module load GLPK/5.0 || true
module load cairo/1.17.8 || true
module load freetype/2.13.0 || true
module load libwebp/1.3.1 || true

# Activate renv (best-effort; master_job.sh will also activate)
PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
if [[ -f "${PROJECT_ROOT}/renv/activate.R" ]]; then
  Rscript -e "source('${PROJECT_ROOT}/renv/activate.R'); renv::load(project='${PROJECT_ROOT}')" >/dev/null 2>&1 || true
fi

# Create logs directory
mkdir -p label_transfer_task/logs

# Run the master job script
bash label_transfer_task/scripts/master_job.sh
