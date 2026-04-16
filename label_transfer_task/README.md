# Label Transfer Task

Comprehensive single-cell classifier benchmarking pipeline for evaluating automated cell-type label transfer methods against reference annotations.

## Overview

This task provides a standardized pipeline to:

1. **Prepare data**: Split preprocessed datasets into reference (30%) and query (70%) sets
2. **Run classifiers**: Execute 14+ different classification methods on query cells using reference annotations
3. **Evaluate**: Measure classifier accuracy and compare performance
4. **Analyze**: Aggregate results and generate visualizations

**Based on:** Consistency_metrics_benchmark/OPS/manual_classifiers/

## Quick Start

### Prerequisites

Ensure you have:
- Completed the `data_processing/` pipeline to generate processed datasets
- R with required packages: `dplyr`, `Seurat`, `SingleR`, `xgboost`, `randomForest`, `e1071`, etc.

### Step 1: Prepare Query/Reference Splits

```bash
cd /path/to/ISC_benchmark_reproducibility
Rscript label_transfer_task/R/00_prepare_splits.R
```

This:
- Loads processed datasets from `data/processed/`
- Filters for datasets with ≥10 samples
- Creates `n_replicates` independent sample-level splits
- Splits each into query (70% samples) and reference (30% samples)
- Saves to `data/processed/label_transfer/rep<k>/<dataset_id>/query.rds` and `reference.rds`

### Step 2: Run All Classifiers

#### Option A: Local Sequential Execution

```bash
cd /path/to/ISC_benchmark_reproducibility
Rscript label_transfer_task/R/01_run_classifiers.R
```

This runs all classifiers on all datasets sequentially.

#### Option B: HPC Batch Submission (Recommended)

```bash
cd /path/to/ISC_benchmark_reproducibility
bash label_transfer_task/scripts/submit_classifiers_hpc.sh
```

Options:
```bash
# Submit jobs for all datasets
bash label_transfer_task/scripts/submit_classifiers_hpc.sh

# Submit job for specific dataset
bash label_transfer_task/scripts/submit_classifiers_hpc.sh JoaI

# Submit jobs for multiple datasets
bash label_transfer_task/scripts/submit_classifiers_hpc.sh JoaI Mitchel BCC
```

This submits independent jobs for each dataset, allowing parallel execution.

## Configuration

All parameters are in `config/label_transfer_parameters.yaml`:

```yaml
seed: 22                    # Reproducibility seed
n_cores: 8                  # Parallel workers
n_replicates: 3             # Independent splits (rep1/rep2/...)

data_preparation:
  min_samples: 10           # Datasets with <10 samples are skipped
  prop_query: 0.7           # % of samples for query (30% reference)

classifiers:
  methods:                  # List of classifiers to benchmark
    - SingleR
    - LogisticRegression
    - XGBoost
    - MLP
    - RandomForest
    - SVM
    - kNN
    - NaiveBayes
    - LDA
    - SeuratTransfer
    - DecisionTree
    - scPred
    - Random                 # Baseline
```

## Directory Structure

```
label_transfer_task/
├── README.md                                # This file
├── config/
│   └── label_transfer_parameters.yaml      # Configuration
├── classifiers/
│   └── classifiers.R                       # Classifier implementations
├── R/
│   ├── 00_utils.R                          # Utilities & paths
│   ├── 00_prepare_splits.R                 # Data preparation script
│   ├── 01_run_classifiers.R                # Classification runner
│   ├── 02_plots_tables.R                   # (Optional)
│   └── 03_consistency.R                    # Consistency + F1 metrics
├── scripts/
│   ├── submit_classifiers_hpc.sh           # HPC submission
│   └── submit_hpc.sh                       # (Legacy)
├── _targets.R                               # Targets workflow (optional)
├── _targets.yaml                            # Targets config (optional)
└── logs/                                    # Job logs (created at runtime)
```

## Classifier Details

### Pure R Implementations (14)

1. **SingleR** - Reference-based annotation using correlation
2. **LogisticRegression** - Log-linear discrimination
3. **XGBoost** - Gradient boosting
4. **MLP** - Multilayer perceptron neural network
5. **RandomForest** - Ensemble decision trees
6. **SVM** - Support vector machine
7. **kNN** - k-Nearest neighbors
8. **NaiveBayes** - Probabilistic Bayesian method
9. **LDA** - Linear discriminant analysis
10. **SeuratTransfer** - Seurat label transfer (MapQuery)
11. **DecisionTree** - Single decision tree
12. **scPred** - Single-cell prediction (scPred)
13. **Random** - Random baseline
14. **Ensemble** - Majority voting from top classifiers

## Data Format

### Input: Reference & Query

**File: `data/processed/label_transfer/rep<k>/<dataset_id>/reference.rds`**

```r
reference <- list(
  counts = <dgCMatrix>,        # Gene × cell sparse matrix
  metadata = <data.frame>      # Columns: sample, cell_type, dataset_id
)
```

**File: `data/processed/label_transfer/rep<k>/<dataset_id>/query.rds`**

```r
query <- list(
  counts = <dgCMatrix>,        # Gene × cell sparse matrix
  metadata = <data.frame>      # Columns: sample, cell_type, dataset_id
)
```

### Output: Predictions

**File: `data/processed/label_transfer/rep<k>/<dataset_id>/predictions.rds`**

Updated query metadata with prediction columns:

```r
query_metadata <- readRDS("predictions.rds")
# Columns: sample, cell_type, dataset_id, pred_<classifier_name>, ...
# pred_SingleR, pred_XGBoost, pred_SVM, ..., pred_Ensemble
```

## Dataset Identification Column Mapping

The pipeline uses dataset-specific identification columns:

```r
get_idents_by_prefix <- function(prefix) {
  switch(prefix,
    "JoaI"       = "cell.type",
    "Mitchel"    = "OriginalAnnotationLevel1",
    "Stephenson" = "OriginalAnnotationLevel2",
    "BCC"        = "annotation",
    "LungAtlas"  = "cell_type",
    "ICBAtlas"   = "cell_type",
    "celltype"   # default fallback
  )
}
```

## Execution Features

### Error Handling

- If a classifier fails, execution continues with remaining classifiers
- Failures are logged but don't interrupt the pipeline
- Invalid predictions (wrong length, all NA) marked with accuracy = NA

### Incremental Execution

- Checks if `predictions.rds` exists before processing
- Skips datasets already completed
- Add new datasets and rerun to process only new data
- Useful for resuming interrupted runs

### Ensemble Voting

- After all classifiers run, ensemble combines predictions
- Uses majority voting across non-Random classifiers
- Typically achieves highest accuracy

## Output Interpretation

For each dataset in predictions.rds:

```r
predictions <- readRDS("data/processed/label_transfer/JoaI/predictions.rds")
head(predictions)

# Columns:
# - sample: Sample ID
# - cell_type: Ground truth annotation
# - dataset_id: Dataset identifier
# - pred_SingleR: SingleR prediction
# - pred_XGBoost: XGBoost prediction
# - ... (one column per classifier)
# - pred_Ensemble: Ensemble prediction

# Calculate per-classifier accuracy:
sapply(predictions[, grep("^pred_", colnames(predictions))], function(col) {
  mean(col == predictions$cell_type, na.rm = TRUE)
})

# Visualize accuracy:
accuracy_df <- data.frame(
  classifier = sub("^pred_", "", grep("^pred_", colnames(predictions), value = TRUE)),
  accuracy = sapply(predictions[, grep("^pred_", colnames(predictions))], 
                    function(col) mean(col == predictions$cell_type, na.rm = TRUE))
)
```

## Workflow Integration (Optional)

To integrate with targets workflow:

```bash
# Load all targets
targets::tar_make()

# Run specific target
targets::tar_make(lt_classifier_results)

# View results
targets::tar_read(lt_classifier_results)
```

## Troubleshooting

### Issue: "This dataset has less than 10 samples"

**Solution:** Increase samples in data_processing if possible, or lower `min_samples` in config.

### Issue: Classifier "FAILED"

**Common causes:**
- Insufficient cells in reference
- Singular matrix (all cells identical)
- Package not installed
- Memory exceeded

Check classifier function in `classifiers/classifiers.R` for verbose error messages.

### Issue: Memory exceeded

**Solution:**
- Reduce `n_cores` in config
- Run datasets in parallel batches
- Check if query/reference have appropriate sample sizes

## Reproducibility Notes

- **Seed**: Controlled by `seed` in config (and offset per replicate to generate different splits)
- **Processors**: Set `n_cores` equal to available cores for reproducible results
- **Software versions**: Lock package versions in renv (if used)

## Citation & References

Based on procedures from:
- **Consistency_metrics_benchmark/OPS/manual_classifiers/**
- Reference implementation for benchmark protocols
- Ensures compatibility with results from other benchmark modules

## Contact & Issues

For issues or questions, refer to:
- Project README: `../README.md`
- Benchmark notes: `../ISC_benchmark/README.md`
└── results/               # (gitignored)
    ├── _targets/          # Targets store (never manually modify)
    ├── raw_results/       # Per-classifier outputs: <dataset>_<clf>_rep<N>.rds
    ├── aggregated/        # Parsed results & summary stats
    └── figures/           # Generated plots
```

## Classifiers

| Classifier | Method | Notes |
|------------|--------|-------|
| RF | Random Forest | 500 trees, default params |
| SVM | Support Vector Machine | RBF kernel, probability estimates |
| KNN | k-Nearest Neighbors | k=5, based on PCA distances |

## Quick Start

### Installation

From the repository root:

```bash
Rscript -e 'renv::restore()'
```

### Run Locally (Sequential)

```bash
cd label_transfer_task
Rscript -e 'targets::tar_make()'
```

This will:
- Load query/reference pairs from `../data/processed/label_transfer/rep<k>/`
- Create grid: datasets × classifiers × replicates
- Run all classifiers (~100-200 runs, depending on participation)
- Aggregate results to `results/aggregated/label_transfer_metrics_aggregated.csv`
- Generate summary plots

Check progress:

```bash
Rscript -e 'targets::tar_progress()'
```

Resume if interrupted:

```bash
Rscript -e 'targets::tar_make()'
```

### Run on HPC (SLURM)

Edit `_targets.yaml` with SLURM resources, then submit via `scripts/submit_hpc.sh`.

## Configuration

All settings are in `config/label_transfer_parameters.yaml`. Edit to configure:

- `seed` — Random seed
- `n_cores` — Number of parallel cores
- `n_replicates` — Number of independent splits (rep1/rep2/...)
- `data_preparation.min_samples` — Datasets with fewer samples are skipped
- `data_preparation.prop_query` — Query split proportion (rest is reference)
- `classifiers.methods` — Which classifiers to evaluate

## Output Files

### In `results/aggregated/`

- `label_transfer_metrics_aggregated.csv` — All classifier runs
  - Columns: `dataset_id`, `classifier`, `replicate`, `accuracy`, `balanced_accuracy`, ...

- `label_transfer_summary_stats.csv` — Summary by classifier

### In `results/figures/`

- Generated plots for paper figures (Figure 5)

## Troubleshooting

### "Label-transfer data directory not found"

Ensure `data_processing/` has been run first to generate query/reference pairs.

### Some datasets missing

Check `../data/processed/<dataset>.rds.yaml` — if `label_transfer: false`, that dataset is excluded.

## Integration with scTypeEval

Consistency metrics are computed via scTypeEval in `R/03_consistency.R`.

## Related

- **Data Processing:** `../data_processing/` — Prepare raw datasets and generate query/reference pairs
- **ISC Benchmark:** `../ISC_benchmark/` — Independent ISC consistency benchmark
- **Main Paper:** Merges ISC + label_transfer results for final figures
