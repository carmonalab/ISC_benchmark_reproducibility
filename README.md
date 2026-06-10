# Reproducibility code for *Evaluating cell type annotations in single-cell omics in the absence of ground truth*

Reproducible workflows used to generate the analyses and figures for the manuscript.

This repository is organized around four main reproducibility layers:
- data preparation from raw single-cell objects
- benchmarking pipelines (inter-sample consistency and label transfer - supervised classification)
- ISC meta-study
- figure notebooks used for manuscript panels

## What You Can Reproduce

- Process all core raw datasets into standardized outputs (`data_processing`)
- Run the ISC benchmark (`ISC_benchmark/`)
- Run saturation analysis on ISC metrics (`ISC_benchmark_saturation/`)
- Run the label-transfer benchmark (supervised classification) (`label_transfer_task/`)
- Run meta-study using top ISC metrics (`Meta_study/`)
- Render manuscript figure notebooks from processed outputs (`Figures_notebooks/`)


## Quick Start

### 1) Environment setup

From repository root:

```bash
Rscript -e 'renv::restore()'
```

### 2) Download raw data

Download raw `.rds` files from Zenodo DOI `10.5281/zenodo.18921437` into:

- `data/raw/`

Dataset-level expected filenames and processing details are documented in:
- `data_processing/README.md`
- `data_processing/config/`

### 3) Run data processing

```bash
cd data_processing
Rscript -e 'targets::tar_make()'
```

Alternative wrappers (local/HPC) are available in `data_processing/scripts/`.

Expected output: processed datasets in `data/processed/`.

## Run Benchmarks

### ISC benchmark

```bash
cd ISC_benchmark
Rscript -e 'targets::tar_make()'
```

Details and HPC submission options:
- `ISC_benchmark/README.md`
- `ISC_benchmark/scripts/`
- `ISC_benchmark/config/`

## Run ISC Saturation Analysis

The saturation pipeline evaluates ranking stability by assessing method rankings across progressively larger subsets of datasets. It measures whether top-performing methods remain robust to dataset selection and predicts ranking convergence at larger dataset numbers.

```bash
cd ISC_benchmark_saturation
Rscript -e 'targets::tar_make()'
```

Main outputs:
- `ISC_benchmark_saturation/results/` — Ranking stability analysis, correlation distributions, and trend curves

Details and optional runner script:
- `ISC_benchmark_saturation/README.md`
- `ISC_benchmark_saturation/scripts/`

### Label-transfer benchmark

The label-transfer pipeline benchmarks supervised classifiers on query/reference splits generated from processed datasets.

Splits are prepared inside the targets pipeline (`lt_prepared_splits` target).

Then run the benchmark:

```bash
cd label_transfer_task
Rscript -e 'targets::tar_make()'
```

Optional: run split generation as a standalone pre-step:

```bash
cd ..
Rscript label_transfer_task/R/00_prepare_splits.R
```

Main aggregated outputs:
- `label_transfer_task/results/aggregated/label_transfer_metrics_aggregated.csv`
- `label_transfer_task/results/aggregated/label_transfer_summary_stats.csv`

Alternative classifier/HPC runners are in:
- `label_transfer_task/R/`
- `label_transfer_task/scripts/`

## Run Meta-study

The meta-study pipeline generates scTypeEval objects from processed datasets and selected annotation columns.

```bash
cd Meta_study
Rscript -e 'targets::tar_make()'
```

Main outputs:
- `Meta_study/output/scTypeEval_objs/`

Details and optional runner script:
- `Meta_study/README.md`
- `Meta_study/scripts/`

## Module Documentation

- `data_processing/README.md`
- `ISC_benchmark/README.md`
- `ISC_benchmark_saturation/README.md`
- `label_transfer_task/README.md`
- `Meta_study/README.md`

## Reproduce Manuscript Figures

Figure notebooks are in `Figures_notebooks/` (for example `1_Figure.Rmd`, `3_Figure.Rmd`, `4_Figure.Rmd`, `5_Figure.Rmd`, `6_Figure.Rmd`).

Typical usage:

```bash
cd Figures_notebooks
Rscript -e 'rmarkdown::render("1_Figure.Rmd")'
```

Notebook respective inputs are listed within notebooks.

## Reproducibility Notes

- Pipeline execution is managed with `targets` for resumable and incremental runs
- Package versions are pinned by `renv.lock`
- Parameters are controlled in YAML files under each module `config/` directory

## Recommended Execution Order

1. `renv::restore()`
2. Download Zenodo raw files to `data/raw/`
3. Run `data_processing/`
4. Run `ISC_benchmark/`
5. Run `ISC_benchmark_saturation/`
6. Run `label_transfer_task/`
7. Run `Meta_study/`
8. Render notebooks in `Figures_notebooks/`

## Troubleshooting

- If a pipeline cannot find paths, run commands from the module directory (`data_processing/`, `ISC_benchmark/`, `ISC_benchmark_saturation/`, `label_transfer_task/`, `Meta_study/`)
- If outputs are partial after interruption, rerun `targets::tar_make()` in the same module to resume
