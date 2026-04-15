#!/bin/bash
#SBATCH --job-name=ISC_data_processing
#SBATCH --partition=private-carmona-gpu
#SBATCH --time=00:45:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=5
#SBATCH --mem=128G
#SBATCH --output=data_processing/logs/processing_%j.log
#SBATCH --error=data_processing/logs/processing_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=josep.garnicacaparros@unige.ch

# ============================================================================
# SLURM Job Submission: ISC Benchmark Dataset Processing via Targets Pipeline
#
# This script submits the data processing job to an HPC cluster using SLURM.
# It runs master_job.sh, which orchestrates the targets-based pipeline.
#
# Before running:
#   1. Edit the SBATCH directives above:
#      - Set YOUR_ACCOUNT to your HPC account
#      - Set YOUR_PARTITION to your compute partition
#      - Set YOUR_EMAIL for notifications
#   2. Ensure raw data is downloaded to data/raw/ (see data_processing/README.md)
#   3. Ensure R has required packages: targets, yaml, dplyr, stringr, Seurat, etc.
#
# Submit the job:
#   sbatch data_processing/scripts/submit_hpc.sh
#
# Monitor:
#   squeue -u $USER
#   tail -f data_processing/logs/processing_<jobid>.log
#
# ============================================================================

# Stop on any error
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
mkdir -p data_processing/logs

# Run the master job script
bash data_processing/scripts/master_job.sh
