# ISC benchmark reproducibility (Garnica et al., 2026)

This repository contains the **reproducible pipelines** used for:

- dataset acquisition + preprocessing
- ISC (Inter-Sample Consistency) benchmarking using `scTypeEval`
- label-transfer benchmarking (manual classifiers)
- figure generation


## Repository layout

- `data_processing/` — download + preprocess datasets from Zenodo
- `ISC_benchmark/` — ISC benchmarking pipeline
- `label_transfer_task/` — label transfer benchmarking pipeline
- `config/` — YAML configuration (datasets + benchmarking parameters)
- `data/` and `results/` — generated artifacts (ignored by git)

## Reproducibility quickstart

### Step 1 — Install dependencies

You need R (>= 4.3 recommended) and the R packages used by the pipeline.

The key dependency is `scTypeEval`.

``` r
# install.packages("remotes")
remotes::install_github("carmonalab/scTypeEval")
```

### Step 2 — Download Zenodo datasets

See `data_processing/README.md` (Zenodo DOI: 10.5281/zenodo.18921437).

Place raw Zenodo `.rds` files into:

- `zenodo/raw/`

### Step 3 — Process datasets

`Rscript data_processing/process_datasets.R`

Outputs:

- `data/processed/isc/*.rds`
- `data/processed/isc/*.rds.yaml` (sidecar metadata, including task participation flags)
- `data/processed/label_transfer/<dataset_id>/{query.rds,reference.rds}`

### Step 4 — Run ISC benchmark (Tasks 1–6)

`Rscript ISC_benchmark/tasks/run_tasks_1_6.R`

Outputs under:

- `results/isc/`

### Step 4b — Run ISC batch/perturbation comparisons (Tasks 7–8)

`Rscript ISC_benchmark/BatchEffect_singleRun/run_tasks_7_8.R`

Outputs under:

- `results/isc/tasks/task7_batch_comparison/output/`
- `results/isc/tasks/task8_perturbation_comparison/output/`

### Step 5 — Run label transfer benchmark

`Rscript label_transfer_task/run_label_transfer_benchmark.R`

Outputs under:

- `results/label_transfer/`

### Step 6 — Reproduce figures

Figure 3 (ISC benchmarking):

`Rscript ISC_benchmark/render_figure3_ISC.R`

Figure 5 (label transfer):

`Rscript -e "rmarkdown::render('label_transfer_task/figure5_label_transfer.Rmd')"`

## Configuration

- `config/datasets.yaml` — maps each dataset reference (as in `specs_datasets.csv`) to Zenodo files + metadata column mappings
- `config/benchmark_parameters.yaml` — all global parameters (seed, cores, downsampling, task parameters)

## Dataset/task mapping

The last four columns of `data_processing/specs_datasets.csv` determine which datasets are used by which task family:

- `ISC-benchmarking` (yes/no): Tasks 1–6
- `Batch comparison` (yes/no): Task 7
- `Perturbation comparison` (yes/no): Task 8
- `Label-Transfer Task` (yes/no): label transfer benchmark
