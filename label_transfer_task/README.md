# Label transfer benchmark

Runs label-transfer benchmarking (manual classifiers) on query/reference splits generated in `data_processing/`.

## Inputs

- `data/processed/label_transfer/<dataset_id>/{query.rds,reference.rds}`

## Output

- `results/label_transfer/` — predictions and summary metrics

## Run

From the repository root:

- `Rscript label_transfer_task/run_label_transfer_benchmark.R`

## Figure 5

From the repository root:

- `Rscript -e "rmarkdown::render('label_transfer_task/figure5_label_transfer.Rmd')"`
