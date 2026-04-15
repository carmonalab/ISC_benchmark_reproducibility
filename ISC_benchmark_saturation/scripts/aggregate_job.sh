#!/bin/bash -l

#SBATCH --account scarmona_til_omics
#SBATCH --mail-type END,FAIL
#SBATCH --mail-user josep.garnicacaparros@unil.ch
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=00:30:00
#SBATCH --export=NONE

set -euo pipefail

# Load required modules (customize for your HPC)
module purge 2>/dev/null || true
module load GCC/12.3.0 2>/dev/null || true
module load R/4.3.2 2>/dev/null || true
module load GLPK/5.0 2>/dev/null || true
module load cairo/1.17.8 2>/dev/null || true
module load freetype/2.13.0 2>/dev/null || true
module load libwebp/1.3.1 2>/dev/null || true

WORKDIR=/scratch/jgarnica/ISC_benchmark_reproducibility/ISC_benchmark_saturation
cd "$WORKDIR"

echo "[HPC] Aggregating cached saturation group results"
Rscript scripts/aggregate_cached_results.R
