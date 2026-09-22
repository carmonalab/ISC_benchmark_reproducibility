#!/bin/bash
#SBATCH --job-name=ISC_label_transfer_between
#SBATCH --partition=public-cpu
#SBATCH --time=94:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --output=label_transfer_task/logs/label_transfer_between_%j.log
#SBATCH --error=label_transfer_task/logs/label_transfer_between_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=josep.garnicacaparros@unige.ch

set -euo pipefail

# Kept in sync with label_transfer_task/scripts/submit_hpc.sh: R 4.5.2 is the only renv
# library tree matching the root renv.lock (e.g. contains anndataR).
module purge
module load GCCcore/10.3.0
module load Python/3.9.5-bare
module load GCC/14.3.0
module load R/4.5.2
module load GLPK/5.0 || true
module load cairo/1.17.8 || true
module load freetype/2.13.0 || true
module load libwebp/1.3.1 || true

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
export RENV_PROJECT="${PROJECT_ROOT}"
export RENV_PROJECT_EXPLICIT="${PROJECT_ROOT}"
export RENV_CONFIG_AUTOLOADER_ENABLED="FALSE"
if [[ -f "${PROJECT_ROOT}/renv/activate.R" ]]; then
  Rscript -e "source('${PROJECT_ROOT}/renv/activate.R'); renv::load(project='${PROJECT_ROOT}')" >/dev/null 2>&1 || true
fi

mkdir -p label_transfer_task/logs

bash label_transfer_task/scripts/master_job_between_datasets.sh
