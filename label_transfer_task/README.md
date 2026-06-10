# Label-Transfer Benchmark Pipeline

Classifier benchmarking pipeline for automated cell-type label transfer reproducibility.

## Scope

This pipeline evaluates multiple classifiers using query/reference splits built from processed datasets.

Main entrypoints:
- `_targets.R`
- `R/00_prepare_splits.R`
- `scripts/master_job.sh`
- `scripts/submit_hpc.sh`

## Inputs

Required:
- Processed datasets from `data_processing/` in `../data/processed/`
- Pipeline parameters in `config/label_transfer_parameters.yaml`

## Quick Start

From repository root:

```bash
Rscript -e 'renv::restore()'
cd label_transfer_task
Rscript -e 'targets::tar_make()'
```

`targets::tar_make()` prepares query/reference splits automatically through the `lt_prepared_splits` target.

Optional standalone split preparation:

```bash
Rscript label_transfer_task/R/00_prepare_splits.R
```

Resume after interruption:

```bash
cd label_transfer_task
Rscript -e 'targets::tar_make()'
```

Optional direct classifier runner:

```bash
cd ..
Rscript label_transfer_task/R/01_run_classifiers.R
```

## HPC Execution

Local/HPC wrapper (runs targets pipeline):

```bash
bash label_transfer_task/scripts/master_job.sh
```

SLURM submission:

```bash
sbatch label_transfer_task/scripts/submit_hpc.sh
```

## Outputs

Aggregated outputs:
- `results/aggregated/label_transfer_metrics_aggregated.csv`
- `results/aggregated/label_transfer_summary_stats.csv`

Figures:
- `results/figures/`

Split inputs and prediction artifacts are written under:
- `../data/processed/label_transfer/rep<k>/<dataset_id>/`

## Configuration

Edit `config/label_transfer_parameters.yaml`:
- `seed`
- `n_cores`
- `n_replicates`
- `data_preparation.min_samples`
- `data_preparation.prop_query`
- `classifiers.methods`

## Troubleshooting

- Missing split files: rerun `targets::tar_make()` (it regenerates splits via `lt_prepared_splits`) or run `R/00_prepare_splits.R` manually.
- Empty/partial results: rerun `targets::tar_make()` to resume.
- Classifier failures: check package availability in the project `renv` environment.

## Related

- Root workflow: `../README.md`
- Data preprocessing: `../data_processing/README.md`
- ISC benchmark: `../ISC_benchmark/README.md`
