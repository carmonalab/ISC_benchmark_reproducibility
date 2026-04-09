#!/bin/bash
#SBATCH --job-name=isc-benchmark
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --partition=normal
#SBATCH --output=slurm-%j.log

# ISC Benchmark pipeline on HPC (SLURM)
#
# Usage:
#   sbatch scripts/submit_hpc.sh
#
# This submits the ISC benchmark pipeline with parallelization via SLURM

cd "$(dirname "$(readlink -f "$0")")/.." || exit 1

# Load R module (adjust for your HPC)
# module load R/4.3

echo "[$(date)] Starting ISC benchmark pipeline on HPC..."

# Run targets with SLURM future plan
Rscript - <<'EOF'
library(targets)
library(future.batchtools)

# Configure SLURM future plan
plan(
  future.batchtools::batchtools_slurm,
  workers = 16,  # Adjust based on available cluster resources
  resources = list(
    n_cores = 4,
    memory_gb = 32,
    walltime = "04:00:00",
    partition = "normal"
  )
)

# Run targets workflow
tar_make()

# Print summary
cat("\n")
tar_summary()
EOF

echo "[$(date)] Pipeline complete."
