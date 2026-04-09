#!/bin/bash
#SBATCH --job-name=ISC_data_processing
#SBATCH --account=<YOUR_ACCOUNT>
#SBATCH --partition=<YOUR_PARTITION>
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --output=data_processing/logs/processing_%j.log
#SBATCH --error=data_processing/logs/processing_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=<YOUR_EMAIL>

# ============================================================================
# SLURM Job Submission: ISC Benchmark Dataset Processing
#
# Before running:
#   1. Edit the SBATCH directives above:
#      - Set YOUR_ACCOUNT to your HPC account
#      - Set YOUR_PARTITION to your compute partition
#      - Set YOUR_EMAIL for notifications
#   2. Ensure raw data is downloaded to data/raw/ (see data_processing/README.md)
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
module load R/4.3  # Adjust version as needed
# module load biology/seurat-4.3  # If available

# Create logs directory
mkdir -p data_processing/logs

# Run the master job script
bash data_processing/scripts/master_job.sh
