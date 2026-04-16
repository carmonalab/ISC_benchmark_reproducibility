#!/bin/bash
# master_job.sh
#
# HPC Entry Point for ISC Benchmark
#
# This script:
# 1. Verifies prerequisites (R, required packages)
# 2. Sets up logging directory
# 3. Submits targets workflow to HPC cluster via SLURM
#
# Usage:
#   cd ISC_benchmark
#   bash master_job.sh
#
# Configuration:
#   - Logi_dir: logs/
#   - Job name: ISC_benchmark_reproducibility
#   - Partition: can be set via SLURM_PARTITION env var
#

set -euo pipefail

# ============================================================================
# CONFIGURATION
# ============================================================================

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
PROJECT_DIR="$( cd "$SCRIPT_DIR/.." && pwd )"
LOG_DIR="$SCRIPT_DIR/logs"
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
LOG_FILE="$LOG_DIR/isc_benchmark_$TIMESTAMP.log"

# HPC Configuration (can override with environment variables)
SLURM_PARTITION="${SLURM_PARTITION:-cpu}"
SLURM_NODES="${SLURM_NODES:-1}"
SLURM_CPUS="${SLURM_CPUS:-8}"
SLURM_MEM="${SLURM_MEM:-64G}"
SLURM_TIME="${SLURM_TIME:-24:00:00}"
MAX_PARALLEL_JOBS="${MAX_PARALLEL_JOBS:-50}"

# ============================================================================
# FUNCTIONS
# ============================================================================

log_message() {
  local level=$1
  shift
  local message="$@"
  echo "[$(date +'%Y-%m-%d %H:%M:%S')] [$level] $message" | tee -a "$LOG_FILE"
}

check_prerequisites() {
  log_message "INFO" "Checking prerequisites..."
  
  # Check if R is available
  if ! command -v R &> /dev/null; then
    log_message "ERROR" "R is not installed or not in PATH"
    exit 1
  fi
  
  # Check R packages
  local r_check=$(R --slave -e "
    packages <- c('targets', 'yaml', 'dplyr', 'tidyr', 'Seurat', 'scTypeEval', 'BiocParallel')
    missing <- packages[!sapply(packages, function(x) requireNamespace(x, quietly=TRUE))]
    if (length(missing) > 0) {
      cat('Missing:', paste(missing, collapse=', '))
      quit('no', 1)
    } else {
      cat('OK')
      quit('no', 0)
    }
  " 2>&1)
  
  if [[ "$r_check" != *"OK"* ]]; then
    log_message "WARNING" "Some R packages may be missing. Continuing anyway..."
  else
    log_message "INFO" "All required R packages found"
  fi
  
  # Check config files
  if [[ ! -f "$SCRIPT_DIR/config/isc_benchmark_parameters.yaml" ]]; then
    log_message "ERROR" "Config file not found: config/isc_benchmark_parameters.yaml"
    exit 1
  fi
  
  if [[ ! -f "$SCRIPT_DIR/config/dataset_idents.yaml" ]]; then
    log_message "ERROR" "Config file not found: config/dataset_idents.yaml"
    exit 1
  fi
  
  # Check processed data
  if [[ ! -d "$PROJECT_DIR/data/processed" ]]; then
    log_message "ERROR" "Processed data directory not found: data/processed"
    log_message "ERROR" "Run data_processing pipeline first"
    exit 1
  fi
  
  log_message "INFO" "All prerequisites satisfied"
}

setup_logging() {
  log_message "INFO" "Setting up logging..."
  
  mkdir -p "$LOG_DIR"
  log_message "INFO" "Log file: $LOG_FILE"
  log_message "INFO" "Working directory: $(pwd)"
  log_message "INFO" "Project directory: $PROJECT_DIR"
}

submit_to_slurm() {
  log_message "INFO" "Submitting to SLURM cluster..."
  log_message "INFO" "Configuration: partition=$SLURM_PARTITION, nodes=$SLURM_NODES, cpus=$SLURM_CPUS, mem=$SLURM_MEM, time=$SLURM_TIME"
  
  # Create SLURM submission script
  local submit_script="$LOG_DIR/submit_$TIMESTAMP.sh"
  
  cat > "$submit_script" << 'SLURM_SCRIPT'
#!/bin/bash
#SBATCH --job-name=ISC_benchmark
#SBATCH --partition=PARTITION
#SBATCH --nodes=NODES
#SBATCH --cpus-per-task=CPUS
#SBATCH --mem=MEM
#SBATCH --time=TIME
#SBATCH --output=LOGDIR/slurm_%j.log
#SBATCH --error=LOGDIR/slurm_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=USER@example.com

set -euo pipefail

# Load modules if needed (uncomment for your system)
# module load R/4.3.0
# module load intel/2023

cd WORKDIR

echo "========== Job Started =========="
echo "Time: $(date)"
echo "Hostname: $(hostname)"
echo "Working directory: $(pwd)"
echo "========== Job Started =========="

# Run targets pipeline
R --slave << 'R_SCRIPT'
library(targets)

# Change to ISC_benchmark directory
setwd("ISC_benchmark")

# Run all targets
message("Starting targets::tar_make() at ", Sys.time())
tar_make()
message("Finished targets::tar_make() at ", Sys.time())

# Cleanup
gc()
invisible(0)
R_SCRIPT

echo "========== Job Finished =========="
echo "Time: $(date)"
echo "========== Job Finished =========="
SLURM_SCRIPT
  
  # Substitute configuration variables
  sed -i "s|PARTITION|$SLURM_PARTITION|g" "$submit_script"
  sed -i "s|NODES|$SLURM_NODES|g" "$submit_script"
  sed -i "s|CPUS|$SLURM_CPUS|g" "$submit_script"
  sed -i "s|MEM|$SLURM_MEM|g" "$submit_script"
  sed -i "s|TIME|$SLURM_TIME|g" "$submit_script"
  sed -i "s|LOGDIR|$LOG_DIR|g" "$submit_script"
  sed -i "s|WORKDIR|$SCRIPT_DIR|g" "$submit_script"
  sed -i "s|USER@example.com|${USER}@$(hostname -f)|g" "$submit_script"
  
  chmod +x "$submit_script"
  
  # Submit to SLURM
  log_message "INFO" "Submitting script: $submit_script"
  local job_id=$(sbatch --parsable "$submit_script")
  
  log_message "INFO" "Job submitted successfully!"
  log_message "INFO" "Job ID: $job_id"
  log_message "INFO" "Check status: squeue -j $job_id"
  log_message "INFO" "Check logs: tail -f $LOG_DIR/slurm_${job_id}.log"
  
  return 0
}

# ============================================================================
# MAIN
# ============================================================================

log_message "INFO" "========== ISC Benchmark HPC Submission =========="
log_message "INFO" "Timestamp: $TIMESTAMP"

setup_logging
check_prerequisites

# Try sbatch, fall back to local execution if not available
if command -v sbatch &> /dev/null; then
  submit_to_slurm
  log_message "INFO" "========== Submission Complete =========="
else
  log_message "WARNING" "sbatch not found. Running locally instead..."
  cd "$SCRIPT_DIR"
  R --slave << 'R_SCRIPT'
library(targets)
message("Starting local targets::tar_make()")
tar_make()
message("Finished!")
invisible(0)
R_SCRIPT
fi
