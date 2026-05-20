# ISC Resource Benchmark

This pipeline benchmarks the runtime in milliseconds, peak memory in MB, and CPU usage of each
ISC method combination defined by:

- one dissimilarity method
- one consistency metric
- one dataset
- one annotation column

It uses the same processed datasets and annotation mappings as the main ISC benchmark,
but stores one output file per ISC combination under `resources/output/` so reruns stay
incremental.

## Run locally

```bash
cd resources
Rscript -e 'targets::tar_make(callr_function = NULL)'
```

Restrict to one dataset:

```bash
cd resources
RESOURCE_DATASET_ID=JoaI_CRC-SG1_Normal Rscript -e 'targets::tar_make(callr_function = NULL)'
```

## Submit to HPC

```bash
sbatch resources/scripts/submit_hpc.sh
```

This submits one job per dataset. Each job uses its own targets store under
`ISC_benchmark_resources/_targets/store_<dataset>` and writes per-method `.rds` files to
`ISC_benchmark_resources/output/<dataset>/<ident>/`.