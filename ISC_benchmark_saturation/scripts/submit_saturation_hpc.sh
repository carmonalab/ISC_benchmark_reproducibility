#!/bin/bash
set -euo pipefail

cd "$(dirname "$0")/.."

N_GROUPS=$(Rscript -e 'cfg <- yaml::read_yaml("config/saturation_parameters.yaml"); cat(length(cfg$analysis$group_sizes))')

echo "Submitting saturation array job with ${N_GROUPS} group sizes"
ARRAY_JOB_ID=$(sbatch --parsable --array=1-"${N_GROUPS}" scripts/saturation_array_job.sh)
echo "Submitted array job: ${ARRAY_JOB_ID}"

AGG_JOB_ID=$(sbatch --parsable --dependency=afterok:"${ARRAY_JOB_ID}" scripts/aggregate_job.sh)
echo "Submitted aggregation job: ${AGG_JOB_ID}"
