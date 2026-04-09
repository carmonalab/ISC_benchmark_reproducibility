# ISC Benchmark Pipeline - Refactored Architecture

## Overview

The ISC (Inter-Sample Consistency) benchmark pipeline has been refactored to follow the proven data_processing pipeline architecture with integrated scTypeEval functions replacing all placeholders.

### Key Improvements

✓ Full scTypeEval integration (no more placeholders)
✓ Targets workflow with dynamic targets per dataset
✓ YAML-based configuration (versionable, reproducible)
✓ 8-task progression from sensitivity to robustness testing
✓ Explicit error handling at each stage
✓ Modular helper functions with clear responsibilities
✓ Batch effect and biological perturbation tests (newly implemented)

---

## Task Progression (8 Tasks)

Tasks are ordered to progress logically from sensitivity tests to robustness evaluations:

### Sensitivity Tests (Tasks 1-2)

1. **missclassify** - Sensitivity to cell type signal degradation
   - Implements: `wr.missclasify()` from scTypeEval
   - Perturbation: Shuffle ~10% cell type labels within samples
   - Type: **Monotonic** (metrics degrade with increasing noise)
   - Rates tested: [1, 0.75, 0.5, 0.25, 0] (proportion of correct labels kept)

2. **SplitCelltype** - Sensitivity to cell type over-partitioning
   - Implements: `wr.splitCellType()` from scTypeEval
   - Perturbation: Artificially split most abundant cell type (50% split rate)
   - Type: **Constant** (metrics should be robust to artificial subtypes)
   - Tests ability to detect when a single cell type is split

### Robustness Tests - Annotation & Data (Tasks 3-6)

3. **Nct** - Robustness to annotation granularity
   - Implements: `wr.Nct()` from scTypeEval
   - Perturbation: Randomly downsample to 50% of cell types
   - Type: **Constant** (tests coarse vs fine annotations)
   - Measures consistency stability with different cell type resolution levels

4. **cellular_complexity** - Robustness to cellular complexity
   - **NEW:** Custom implementation (not in original scTypeEval)
   - Perturbation: Stratify cells by gene expression variance (high vs low complexity)
   - Type: **Constant** (tests robustness across tissue complexity tiers)
   - Compares consistency in low-variability vs high-variability populations

5. **Nsamples** - Robustness to dataset size (samples)
   - Implements: `wr.NSamples()` from scTypeEval
   - Perturbation: Randomly downsample to 50% of samples
   - Type: **Constant** (tests missing data robustness)
   - Rates tested: [1, 0.9, 0.7, 0.5]

6. **NCell** - Robustness to dataset size (cells per cell type)
   - Implements: `wr.NCell()` from scTypeEval
   - Perturbation: Reduce cells per sample (shallow sequencing)
   - Type: **Constant** (tests sequencing depth robustness)
   - Rates tested: [1, 0.75, 0.5, 0.25]

### Robustness Tests - Experimental Factors (Tasks 7-8)

7. **batch_effects** - Robustness to batch effects
   - **NEW:** Custom implementation
   - Perturbation: Compare within-batch vs cross-batch consistency
   - Type: **Constant** (tests batch correction robustness)
   - Requires: Batch/site/dataset column in metadata
   - Output: Within-batch and cross-batch metrics

8. **biological_perturbations** - Robustness to biological perturbations
   - **NEW:** Custom implementation
   - Perturbation: Compute consistency per biological condition
   - Type: **Constant** (tests biological signal robustness)
   - Requires: Condition column in metadata
   - Evaluates ISC across treatment/control conditions

---

## Pipeline Architecture

Follows the same structure as `data_processing/`:

```
ISC_benchmark/
├── _targets.R                                  # Targets workflow (dynamic targets)
├── config/
│   └── isc_benchmark_parameters.yaml          # Configuration (8 tasks, metrics, thresholds)
├── R/
│   ├── 01_isc_tasks.R                        # Task definitions & orchestration
│   └── isc_benchmark_helpers.R               # Helper functions wrapping scTypeEval
└── README.md                                   # This file
```

---

## Configuration (YAML)

File: `config/isc_benchmark_parameters.yaml`

### Structure

```yaml
seed: 22                                       # Reproducibility seed
n_cores: 8                                    # Parallel processing

run:                                          # Enable/disable each task
  missclassify: true
  SplitCelltype: true
  Nct: true
  cellular_complexity: true
  Nsamples: true
  NCell: true
  batch_effects: true
  biological_perturbations: true

common:                                       # Parameters shared by all tasks
  normalization_method: "Log1p"              # Log1p or Seurat normalization
  reduction: true                             # Apply PCA reduction
  n_dims: 30                                  # Number of PCA dimensions
  
  dissimilarity_methods:                     # Methods for computing consistency
    - "WasserStein"
    - "Pseudobulk:Euclidean"
    - "Pseudobulk:Cosine"
    - "Pseudobulk:Pearson"
    - "BestHit:Match"
    - "BestHit:Score"
  
  intrinsic_validation_metrics:              # Quality metrics
    - "silhouette"
    - "2label.silhouette"
    - "NeighborhoodPurity"
    - "ward.PropMatch"
    - "Orbital.medoid"
    - "Average.similarity"

task_missclassify:                           # Task 1 parameters
  rates: [1, 0.75, 0.5, 0.25, 0]
  replicates: 3
  shuffle_rate: 0.1

task_SplitCelltype:                          # Task 2 parameters
  rates: [1, 0.5]
  replicates: 3
  split_rate: 0.5

# ... (6 more task sections with task-specific parameters)
```

---

## Key Functions

### isc_benchmark_helpers.R

**Data Preparation:**
- `prepare_scTypeEval_object(obj, ident_col, config)` - Format Seurat for scTypeEval

**Task Runners (wrap scTypeEval wr_* functions):**
- `run_task_missclassify()` → `scTypeEval:::wr.missclasify()`
- `run_task_SplitCelltype()` → `scTypeEval:::wr.splitCellType()`
- `run_task_Nct()` → `scTypeEval:::wr.Nct()`
- `run_task_cellular_complexity()` → Custom (variance-based stratification)
- `run_task_Nsamples()` → `scTypeEval:::wr.NSamples()`
- `run_task_NCell()` → `scTypeEval:::wr.NCell()`
- `run_task_batch_effects()` → Custom (within/cross-batch comparison)
- `run_task_biological_perturbations()` → Custom (per-condition analysis)

**Metric Processing:**
- `extract_task_metrics(wr_result, task_name, metric_type)` - Parse wr_* outputs
- `save_task_results(...)` - Persist results and metadata

### 01_isc_tasks.R

**Main Orchestrators:**
- `run_isc_benchmark_on_dataset(dataset_id, ident_col, task_name, ...)` 
  - Execute single task on single dataset
  - Loads → prepares → runs task → extracts metrics → saves results
  
- `run_isc_benchmark_full(dataset_id, dataset_path, ident_cols, ...)`
  - Execute all active tasks on all cell type annotations
  - Aggregates results into unified data frame
  
- `aggregate_isc_results(result_files, output_file)`
  - Combine results across datasets into summary statistics

---

## Execution Modes

### 1. Targets Workflow (Recommended)

```bash
cd ISC_benchmark
targets::tar_make()                          # Run all pending targets
targets::tar_status()                        # Check execution status
targets::tar_read(summary_all)              # Load aggregated results
```

**Why Targets?**
- Incremental execution (skip completed tasks)
- Automatic dependency tracking
- Integrated caching
- Parallel-ready configuration

### 2. Standalone Script (Not Yet Implemented)

```bash
# Following data_processing pattern:
Rscript run_isc_benchmark.R --dry-run
Rscript run_isc_benchmark.R
```

---

## Output Structure

Results saved to: `results/isc_benchmark/`

```
results/isc_benchmark/
├── all_results.rds                          # Full aggregated results (data frame)
├── dataset_id_1/
│   ├── {task}_{ident}/
│   │   ├── dataset_task_ident_metrics.rds   # Metrics data frame
│   │   ├── dataset_task_ident_metadata.yaml # Config snapshot + seed + timestamp
│   │   └── dataset_task_ident_wrobj.rds     # Full wr_* object (optional)
│   └── ... (other task_ident combinations)
└── dataset_id_2/
    └── ... (same structure)
```

### Metrics Files

Data frames containing:
- `task` - Task identifier
- `ident` - Cell type annotation column name
- `dataset_id` - Dataset identifier
- `status` - "success" or "failed"
- `error` - Error message (if failed)
- Task-specific metric columns (varies)

### Metadata Files (YAML)

Contains:
- Task name, dataset, cell type annotation
- Seed, timestamp, n_cores
- Full configuration snapshot
- R session info for reproducibility

---

## Changes from Original Implementation

| Aspect | Old (Placeholder) | New (Refactored) |
|--------|------------------|------------------|
| **wr_* functions** | Not called | Fully integrated |
| **Tasks 4, 7, 8** | Not implemented | Fully implemented |
| **Configuration** | Hard-coded in R | YAML-based |
| **Workflow** | No targets | Full targets pipeline |
| **Error handling** | Basic try-catch | Explicit per-stage handling |
| **Reproducibility** | Unclear seed propagation | Clear seed management |
| **Modularity** | Mixed concerns | Separated: helpers vs orchestration |

---

## Reference Files & Attribution

Original implementations used as templates:
- `data_processing/_targets.R` - Targets workflow pattern
- `data_processing/config/processing_parameters.yaml` - YAML structure
- `data_processing/R/data_processing_helpers.R` - Helper organization
- `Consistency_metrics_benchmark/tasks/Metrics_Benchmarking.R` - Task implementations
- `scTypeEval/inst/benchmarking/Metrics_benchmarking.R` - wr_* function definitions

---

## Roadmap

- [ ] Create standalone script `run_isc_benchmark.R`
- [ ] Add visualization functions (task-specific plots)
- [ ] Extend `cellular_complexity` with true SCAP stratification
- [ ] Implement `batch_correction` task (ComBat, Harmony, etc.)
- [ ] Add `mergeCT` compatibility test
- [ ] Cross-dataset comparison aggregation
- [ ] Sensitivity analysis on perturbation rates
- [ ] Benchmark execution guide documentation

---

## Troubleshooting

### "No processed datasets found"
- Run `../data_processing/_targets.R` first to generate preprocessed files
- Check: `ls ../data/processed/isc/`

### Tasks failing silently
- Check result `.yaml` files: `cat results/isc_benchmark/dataset_id/task_ident/metadata.yaml`
- Review error messages in `.rds` result files

### Memory issues with large datasets
- Reduce `n_cores` in config (parallel overhead)
- Set `save_wr_objects: false` to skip saving large wr_* objects

---

## Contact & Support

For questions about:
- **Task definitions**: See task-specific documentation in `isc_benchmark_parameters.yaml`
- **scTypeEval integration**: Consult scTypeEval package documentation
- **Pipeline architecture**: See `../data_processing/README.md` for reference pattern

Each task is implemented as a function in the `ISC_TASK_CATALOG` list in `R/01_isc_tasks.R`. The implementation follows the logic of scTypeEval's corresponding `wr_*` wrapper functions, but adapted for Seurat objects:

#### Helper Functions (from scTypeEval/inst/benchmarking/)

- `rand_shuffling_group(vector, group, rate=0.1, seed)` — Shuffles a proportion of labels within groups (Task 1, 6)
  - Source: `Metrics_benchmarking.R:22`
- `rand_split_group(vector, group, rate=0.5, celltype, seed)` — Splits a cell type across samples at given rate (Task 6)
  - Source: `Metrics_benchmarking.R:53`
- `downsample_factor_level(df, factor_col, level, threshold, seed)` — Downsamples a factor level to threshold (utility)
  - Source: `Metrics_benchmarking.R:1`

#### Task Runtime Pattern

Each task execution (in `run_isc_benchmark_on_dataset()`):
1. Loads a processed Seurat object from `../data/processed/isc/<dataset>.rds`
2. Applies the perturbation (e.g., shuffle labels) to get `obj_perturbed`
3. Calls `compute_isc_metrics(obj_original, obj_perturbed, ...)`
4. Saves results + metadata + session info to `results/task_results/`

#### Metric Computation (From wr_* Functions)

The `compute_isc_metrics()` function follows the execution pattern from all `wr_*` functions:

```r
# Pattern from wr_nsamples, wr_nct, etc (Metrics_benchmarking.R)
wrapper_dissimilarity(
  sc,
  ident = ident,
  sample = sample,
  dissimilarity_method = c("WasserStein", "Pseudobulk:Euclidean", ...),
  ...
) → get_consistency(sc)  # Returns metrics dataframe
```

Actual call sequence:
1. `run_processing_data()` — Normalize, filter by min_samples/min_cells
2. `add_gene_list()` — Select gene markers (HVG or custom)
3. `run_pca()` — Dimensionality reduction (optional)
4. `run_dissimilarity()` — Compute dissimilarities (WasserStein, Pseudobulk variants, reciprocal classification)
5. `get_consistency()` — Extract consistency scores (silhouette, neighborhood purity, etc.)

### Code Origin & Attribution

**Primary source:**
- `scTypeEval/inst/benchmarking/Metrics_benchmarking.R` — Complete `wr_*` implementations
  - `wr_missclasify` (line 227)
  - `wr_nsamples` (line 377)
  - `wr_nct` (line 532)
  - `wr_ncell` (line 678)
  - `wr_merge_ct` (line 849)
  - `wr_split_cell_type` (line 1018)
  - `wrapper_dissimilarity` (line 75)

**Helper functions:**
- `scTypeEval/inst/benchmarking/Metrics_benchmarking.R` (lines 1–120)
  - `rand_shuffling_group`, `rand_split_group`, `downsample_factor_level`, etc.
- `scTypeEval/inst/benchmarking/assays_utils.R` (utility functions)

**Legacy implementations:**
- `Consistency_metrics_benchmark/tasks/Metrics_Benchmarking.R` (original test code)
- `Consistency_metrics_benchmark/BatchEffect_singleRun/BatchEffect_singleRun.R` (batch tasks)

**Note:** All function names have been adapted to snake_case following R conventions (e.g., `rand_shuffling_group` instead of camelCase variants).

**Note:** Tasks 7–8 (batch effect) can be enabled in `config/isc_parameters.yaml` for future expansion.

## Quick Start

### Installation

From the repository root:

```bash
# Install R dependencies (if needed)
Rscript -e 'renv::restore()'

# Or manually:
Rscript -e 'install.packages("targets"); install.packages("tarchetypes")'
```

### Run Locally (Sequential)

```bash
cd ISC_benchmark
Rscript -e 'targets::tar_make()'
```

This will:
- Load processed datasets from `../data/processed/isc/`
- Create grid: datasets × tasks × replicates
- Run all ~400 ISC benchmark tasks (parallelizable; see below)
- Aggregate results to `results/aggregated/isc_metrics_aggregated.csv`
- Generate summary figures

Check progress:

```bash
Rscript -e 'targets::tar_progress()'
```

Resume if interrupted:

```bash
Rscript -e 'targets::tar_make()'  # Picks up where it left off
```

### Run on HPC (SLURM)

Edit `_targets.yaml` with SLURM resource requirements:

```yaml
controller: future
future:
  plan: future.batchtools::batchtools_slurm
  workers: 16              # Max parallel jobs
  
resources:
  slurm_time: "04:00:00"
  slurm_cpus_per_task: 4
  slurm_memory_gb: 32
  slurm_partition: "normal"
```

Then submit via `scripts/submit_hpc.sh`.

## Configuration

All settings are in `config/isc_parameters.yaml`. Edit to configure:

- `seed` — Random seed (for reproducibility)
- `ncores` — Number of parallel cores
- `isc_tasks.run` — Enable/disable individual tasks (1-6)
- `isc_tasks.common` — Metric computation settings (reduction, ndim, metrics, dissimilarity methods)
- `grid.n_replicates` — Number of replicates per task-dataset pair
- `processing.min_samples` / `processing.min_cells` — Filters for consistency metrics

## Output Files

### In `results/aggregated/`

- `isc_metrics_aggregated.csv` — Combined results from all runs
- `isc_summary_stats.csv` — Summary statistics

### In `results/figures/`

- Generated plots for paper figures

## Troubleshooting

### "Could not find project root"

Ensure you run from `ISC_benchmark/` or that the project root contains `.Rproj` and `renv.lock` files.

### Some tasks skipped

Check `../data/processed/isc/<dataset>.rds.yaml` — if `isc: false`, that dataset is excluded.

## Integration with scTypeEval

This pipeline uses metrics and utilities from the `scTypeEval` package:

```r
remotes::install_github("carmonalab/scTypeEval")
```

### Metric Computation Workflow

The full ISC benchmark workflow (placeholder for implementation):

1. **Task perturbation** (`R/01_isc_tasks.R`)
   - Load processed dataset
   - Apply task perturbation (e.g., misclassify 10% of cells)
   - Keep original for comparison

2. **Consistency metric computation** (`R/02_scTypeEval_helpers.R`)
   - `run_isc_dissimilarity()` — Compute pairwise dissimilarities
   - `compute_isc_metrics()` — Extract consistency scores
   - Supported metrics: silhouette, neighborhood_purity, average_similarity, etc.

3. **Quantify degradation**
   - Compare original vs perturbed consistency
   - Compute drop: `original_consistency - perturbed_consistency`
   - Evaluate robustness of metrics to perturbation

### Implemented Helper Functions

From `scTypeEval/inst/benchmarking/`:

- `rand_shuffling_group()` — Shuffle labels within groups (10% rate)
- `rand_split_group()` — Split a cell type across samples (50% rate)
- `downsample_factor_level()` — Downsample to threshold
- `monotonicity_score()` — Evaluate monotonic degradation
- `consistency_drop()` — Quantify split effects
- `fit_constant()` — Fit constant model to metric trajectory

### Next Steps: Full Integration

To compute actual metrics, implement in `R/02_scTypeEval_helpers.R`:

1. Call `scTypeEval::run_processing_data()` on both original/perturbed objects
2. Call `scTypeEval::run_pca()` for dimensionality reduction
3. Call `scTypeEval::run_dissimilarity()` with specified methods (Wasserstein, Pseudobulk, etc.)
4. Call `scTypeEval::get_consistency()` to extract metric values
5. Aggregate into per-task result tables

Example (pseudocode):

```r
# Compute dissimilarity metrics
sc_original <- run_isc_dissimilarity(obj_original, ...)
sc_perturbed <- run_isc_dissimilarity(obj_perturbed, ...)

# Extract consistency scores
meta_original <- extract_isc_consistency_results(sc_original, metric = "silhouette")
meta_perturbed <- extract_isc_consistency_results(sc_perturbed, metric = "silhouette")

# Quantify degradation
drop <- quantify_metric_degradation(meta_original, meta_perturbed)
```

## Related

- **Data Processing:** `../data_processing/` — Prepare raw Zenodo datasets
- **Label Transfer:** `../label_transfer_task/` — Independent label-transfer benchmark
- **Paper:** Merges ISC + label_transfer results for final figures
- **scTypeEval:** [GitHub link to scTypeEval repo]
