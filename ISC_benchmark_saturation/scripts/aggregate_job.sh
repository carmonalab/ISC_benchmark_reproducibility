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

module load r-light/4.3.3

WORKDIR=/scratch/jgarnica/ISC_benchmark_reproducibility/ISC_benchmark_saturation
cd "$WORKDIR"

echo "[HPC] Aggregating cached saturation group results"
Rscript scripts/aggregate_cached_results.R
