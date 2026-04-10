#!/bin/bash
set -euo pipefail

cd "$(dirname "$0")/.."

GROUP_SIZES=$(Rscript -e 'cfg <- yaml::read_yaml("config/saturation_parameters.yaml"); cat(paste(cfg$analysis$group_sizes, collapse=" "))')

for k in $GROUP_SIZES; do
  echo "[LOCAL] Running group size $k"
  Rscript scripts/run_group_size.R "$k"
done

echo "[LOCAL] Aggregating cached results"
Rscript scripts/aggregate_cached_results.R

echo "[LOCAL] Done"
