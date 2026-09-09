# ISC Resource Benchmark

This pipeline benchmarks the runtime in milliseconds, peak memory in MB, and CPU usage of every
tool used by `ISC_benchmark`, ncores = 1, median of 3 replicate runs each:

- every internal scTypeEval `dissimilarity_method` x `int_val_metric` pair
- every enabled external tool (SCCAF, anticor_features, sc-SHC), R or Python

for each dataset / annotation column combination.

Metrics are measured with `/usr/bin/time` wrapping the exact top-level process being timed
(an `Rscript` child for internal/sc-SHC methods, or the `python` interpreter directly for
SCCAF/anticor_features), with BLAS/OMP thread pools pinned to `ncores = 1` via environment
variables. This makes duration, CPU usage and peak memory directly comparable between R and
Python tools.

It uses the same processed datasets and annotation mappings as the main ISC benchmark,
but stores one output file per tool under `output/<dataset>/<ident>/` so reruns stay
incremental.

## Configuration

`config/resource_parameters.yaml` mirrors `ISC_benchmark/config/isc_benchmark_parameters.yaml`
(`common.dissimilarity_method`, `common.int_val_metric`, `external_methods`) so every method
actually exercised by `ISC_benchmark` gets profiled here too.

## Run locally

```bash
cd ISC_benchmark_resources
Rscript -e 'targets::tar_make(callr_function = NULL)'
```

Restrict to one dataset:

```bash
cd ISC_benchmark_resources
RESOURCE_DATASET_ID=JoaI_CRC-SG1_Normal Rscript -e 'targets::tar_make(callr_function = NULL)'
```

## Submit to HPC

```bash
bash scripts/submit_hpc.sh
```

This submits one job per dataset. Each job uses its own targets store under
`ISC_benchmark_resources/_targets/store_<dataset>` and writes per-tool `.rds` files to
`ISC_benchmark_resources/output/<dataset>/<ident>/`. Once every per-dataset job completes,
`submit_hpc.sh` automatically submits `scripts/submit_aggregate.sh` (via
`--dependency=afterok:...`) to append all datasets' outputs into a single combined table at
`results/aggregated_benchmarks.rds` / `.csv`.