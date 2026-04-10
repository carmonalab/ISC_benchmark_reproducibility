#!/bin/bash
# ==============================================================================
# HPC Job Submission Script for Label Transfer Classification
# ==============================================================================
# Usage: ./scripts/submit_classifiers_hpc.sh
#        ./scripts/submit_classifiers_hpc.sh <dataset_id>
#
# Submits individual jobs for each dataset in the label_transfer_task pipeline
# Each job runs all classifiers for one dataset independently
# ==============================================================================

set -e

# Project root
PROJECT_ROOT="/Users/garnica/Documents/Projects/ISC_benchmark_reproducibility"
SCRIPT_DIR="${PROJECT_ROOT}/label_transfer_task"
DATA_DIR="${PROJECT_ROOT}/data/processed/label_transfer"
LOG_DIR="${SCRIPT_DIR}/logs"

# Create log directory
mkdir -p "$LOG_DIR"

# ============================================================================
# FUNCTION: Submit individual dataset job
# ============================================================================

submit_dataset_job() {
    local dataset_id=$1
    local dataset_dir="${DATA_DIR}/${dataset_id}"
    local log_file="${LOG_DIR}/${dataset_id}_classify.log"
    
    if [ ! -d "$dataset_dir" ]; then
        echo "⚠ Skipping $dataset_id: directory not found"
        return 1
    fi
    
    # Check if already completed
    if [ -f "${dataset_dir}/predictions.rds" ]; then
        echo "✓ Skipping $dataset_id: already completed"
        return 0
    fi
    
    echo "→ Submitting job for $dataset_id..."
    
    # Create inline R script
    cat > /tmp/classify_${dataset_id}.R << 'EOF'
#!/usr/bin/env Rscript
suppressPackageStartupMessages({
    library(dplyr)
    library(yaml)
})

# Source utilities
project_root <- Sys.getenv("PROJECT_ROOT", "/Users/garnica/Documents/Projects/ISC_benchmark_reproducibility")
source(file.path(project_root, "R/cli_utils.R"))
source(file.path(project_root, "R/shared_helpers.R"))
source(file.path(project_root, "label_transfer_task/R/00_utils.R"))

# Load classification function
source(file.path(project_root, "label_transfer_task/R/01_run_classifiers.R"))

# Run classification for this dataset
dataset_id <- Sys.getenv("DATASET_ID")
base_dir <- Sys.getenv("DATA_DIR")
params <- load_lt_params()

dataset_dir <- file.path(base_dir, dataset_id)
classify_dataset(dataset_dir, dataset_id, params)
EOF

    # Submit job to HPC (adjust to your HPC system: slurm, pbs, lsf, etc.)
    # This example uses bash/local execution; replace with your scheduler
    
    (
        export PROJECT_ROOT="$PROJECT_ROOT"
        export DATASET_ID="$dataset_id"
        export DATA_DIR="$DATA_DIR"
        
        cd "$PROJECT_ROOT"
        Rscript /tmp/classify_${dataset_id}.R >> "$log_file" 2>&1
    ) &
    
    echo "  → Submitted with log: $log_file (PID: $!)"
}

# ============================================================================
# MAIN EXECUTION
# ============================================================================

if [ ! -d "$DATA_DIR" ]; then
    echo "ERROR: Data directory not found: $DATA_DIR"
    exit 1
fi

echo "=============================================================="
echo "LABEL TRANSFER CLASSIFICATION - HPC JOB SUBMISSION"
echo "=============================================================="
echo "Project: $PROJECT_ROOT"
echo "Data dir: $DATA_DIR"
echo "Log dir: $LOG_DIR"
echo ""

# Get list of datasets
if [ $# -eq 0 ]; then
    # Submit jobs for all datasets
    echo "Submitting jobs for all datasets..."
    for dataset_dir in "$DATA_DIR"/*; do
        if [ -d "$dataset_dir" ] && [ "$(basename "$dataset_dir")" != ".DS_Store" ]; then
            dataset_id=$(basename "$dataset_dir")
            submit_dataset_job "$dataset_id"
        fi
    done
    
    echo ""
    echo "✓ All jobs submitted!"
    echo "Monitor logs in: $LOG_DIR"
else
    # Submit job for specific dataset(s)
    for dataset_id in "$@"; do
        submit_dataset_job "$dataset_id"
    done
fi

# Wait for all background jobs to complete (optional)
# Uncomment to wait for all jobs:
# wait

echo "Done."
