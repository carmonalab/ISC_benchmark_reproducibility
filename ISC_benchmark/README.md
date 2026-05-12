# ISC Benchmark Pipeline

## Overview

This pipeline benchmarks inter-sample consistency (ISC) robustness across 8 tasks using processed single-cell datasets.

It uses:
- `targets` for orchestration and resumability
- `scTypeEval` for consistency computation
- SLURM submission wrappers for HPC execution

Main entrypoints:
- `ISC_benchmark/_targets.R`
- `ISC_benchmark/scripts/submit_hpc.sh`

## What The Pipeline Runs

Tasks are grouped as:

1. Stem-level tasks (single dataset files):
- `missclassify`
- `SplitCelltype`
- `Nct`
- `cellular_complexity`
- `Nsamples`
- `NCell`

2. Family-level tasks (merged family datasets):
- `batch_effects`
- `biological_perturbations`

Behavior for tasks 7-8:
- Pair candidates are derived from `data_processing/config/specs_datasets.csv`.
- Only resolvable pairs for active merged stems are kept.
- For each valid pair:
  - load single-object baseline consistency when available
  - merge the two single objects
  - compute merged consistency
  - return raw consistency rows for single + merged (not ratio/degradation collapse)

## Current Output Semantics

Primary aggregate output:
- `ISC_benchmark/results/all_results.csv`

Important details:
- This is built from in-memory `task_result` returns in `_targets.R`.
- It is cumulative across runs (existing CSV is read, new rows appended, exact duplicates removed).
- It is not a strict filesystem concatenation of `*_metrics.rds` files.

Per-task output files:
- For persisted tasks, metrics are saved as:
  - `results/<dataset_or_family>/<task>_<ident>/<dataset>_<task>_<ident>_metrics.rds`
  - `results/<dataset_or_family>/<task>_<ident>/<dataset>_<task>_<ident>_metadata.yaml`

No-pair handling for tasks 7-8:
- If no resolvable pairs are found, task returns `status=success` with `n_results=0`.
- Persistence is skipped, so no task metrics file/directory is created for that case.

Cache behavior for tasks 7-8:
- Single subset objects may be cached on disk as `subset_sc_*`.
- Merged pair objects (`merged_sc_*`) are computed in memory and are not saved.

## Configuration

Main config:
- `ISC_benchmark/config/isc_benchmark_parameters.yaml`

Key sections:
- `run`: enable/disable each task
- `common`: shared scTypeEval parameters (`verbose`, methods, filtering, etc.)
- `task_*`: task-specific parameters

`verbose` flag:
- Controlled by `common.verbose` (currently `false` unless you change it).

## How To Launch

### 1. Run all tasks for all datasets/families (recommended default)

From `ISC_benchmark/`:

```bash
unset ISC_DATASET_IDS ISC_DATASET_ID ISC_TEST_DATASET ISC_DATASET_FAMILIES ISC_DATASET_FAMILY ISC_TASKS
bash scripts/submit_hpc.sh
```

This submits:
- one stem job per dataset for tasks 1-6
- one family job per family for tasks 7-8

### 2. Run only a subset of stem datasets (tasks 1-6)

```bash
ISC_DATASET_IDS="JoaI_CRC-SG1_Normal,Stephenson_Cambridge_Covid" bash scripts/submit_hpc.sh
```

or

```bash
ISC_DATASET_ID="JoaI_CRC-SG1_Normal" bash scripts/submit_hpc.sh
```

### 3. Run only selected families (tasks 7-8)

```bash
unset ISC_DATASET_IDS ISC_DATASET_ID ISC_TEST_DATASET ISC_TASKS
export ISC_DATASET_FAMILIES="BCC,ICBAtlas,JoaI,LungAtlas,Stephenson"
bash scripts/submit_hpc.sh
```

