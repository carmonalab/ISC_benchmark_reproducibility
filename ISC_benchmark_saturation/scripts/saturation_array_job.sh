#!/bin/bash -l

#SBATCH --account scarmona_til_omics
#SBATCH --mail-type END,FAIL
#SBATCH --mail-user josep.garnicacaparros@unil.ch
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=120G
#SBATCH --time=02:00:00
#SBATCH --export=NONE

set -euo pipefail

module load r-light/4.3.3

WORKDIR=/scratch/jgarnica/ISC_benchmark_reproducibility/ISC_benchmark_saturation
cd "$WORKDIR"

if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
  echo "SLURM_ARRAY_TASK_ID is not set"
  exit 1
fi

GROUP_SIZE=$(Rscript -e 'cfg <- yaml::read_yaml("config/saturation_parameters.yaml"); gs <- cfg$analysis$group_sizes; idx <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID")); if (idx < 1 || idx > length(gs)) stop("Array index out of range"); cat(gs[[idx]])')

echo "[HPC] Running group size ${GROUP_SIZE} (array index ${SLURM_ARRAY_TASK_ID})"
Rscript scripts/run_group_size.R "$GROUP_SIZE"
