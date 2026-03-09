# ISC benchmark

Runs the ISC (Inter-Sample Consistency) benchmarking tasks using the `scTypeEval` package.

Manuscript mapping:

- **Tasks 1–6**: `ISC_benchmark/tasks/run_tasks_1_6.R`
- **Tasks 7–8**: `ISC_benchmark/BatchEffect_singleRun/run_tasks_7_8.R`

The dataset participation for each task family is controlled by the last four columns of `data_processing/specs_datasets.csv`:

- `ISC-benchmarking` → Tasks 1–6
- `Batch comparison` → Task 7
- `Perturbation comparison` → Task 8

## Inputs

- `data/processed/isc/*.rds`
- `config/benchmark_parameters.yaml`

## Output

- `results/isc/tasks/` — task outputs
- `results/isc/resources/` — resource profiling outputs

## Run

From the repository root:

- Tasks 1–6: `Rscript ISC_benchmark/tasks/run_tasks_1_6.R`
- Tasks 7–8: `Rscript ISC_benchmark/BatchEffect_singleRun/run_tasks_7_8.R`

## Figure 3

From the repository root:

- `Rscript ISC_benchmark/render_figure3_ISC.R`
