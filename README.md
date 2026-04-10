# ISC Benchmark Reproducibility (Garnica et al., 2026)

A production-ready, HPC-compatible benchmarking pipeline in R for reproducible single-cell analysis evaluation.

This repository contains **two independent targets-based pipelines**:

1. **ISC Benchmark** — Inter-sample consistency metrics across 6 metadata/count perturbations
2. **Label-Transfer Benchmark** — Manual cell-type classifier evaluation

Both pipelines share preprocessed data but operate independently, enabling parallel execution and clear separation of concerns.

---

## 🚀 Quick Start (Running the Full Pipeline)

### Prerequisites
- R ≥ 4.0 with packages: `targets`, `yaml`, `Seurat`, `dplyr`
- Raw datasets downloaded from Zenodo (see below)

### Step 1: Download Raw Data
Download all `.rds` files from **Zenodo DOI: 10.5281/zenodo.18921437** to `data/raw/`

```bash
cd /path/to/ISC_benchmark_reproducibility
# Or manually download from https://zenodo.org/record/18921437
ls -lh data/raw/*.rds  # Verify downloads
```

### Step 2: Process All Datasets
Submit the data processing job (requires raw data in place):

```bash
# Local processing (serial)
bash data_processing/scripts/master_job.sh

# OR submit to HPC (SLURM)
# First edit data_processing/scripts/submit_hpc.sh to configure your HPC account
sbatch data_processing/scripts/submit_hpc.sh
```

**Once processing completes**, you'll have clean data in `data/processed/isc/` and `data/processed/label_transfer/`

### Step 3: Prepare Label Transfer Data (Optional)
If running the label transfer benchmark, prepare query/reference splits:

```bash
# Prepare label transfer query/reference splits (from project root)
Rscript label_transfer_task/R/00_prepare_splits.R

# This creates data/processed/label_transfer/ with query.rds and reference.rds for each dataset
```

### Step 4: Run Benchmarks
After successful data processing:

```bash
# ISC Benchmark (inter-sample consistency)
cd ISC_benchmark
Rscript -e 'targets::tar_make()'
# Results → ISC_benchmark/results/

# Label Transfer Benchmark (classification, requires Step 3)
cd ../label_transfer_task
Rscript -e 'targets::tar_make()'
# Results → label_transfer_task/results/
```

---

## Repository Structure

```
ISC_benchmark_reproducibility/
├── README.md                          # This file
├── renv.lock                          # Reproducible R environment snapshot
├── ISC_benchmark.Rproj                # RStudio project
│
├── utils/                             # ⭐ Shared utilities
│   ├── shared_helpers.R               # Common paths, config loading, reproducibility
│   ├── project_paths.R                # Path helpers (legacy, keep for now)
│   ├── cli_utils.R                    # Message utilities
│   └── data_processing_helpers.R      # Data I/O helpers
│
├── data/                              # ⭐ Shared processed data (gitignored)
│   ├── raw/                           # Zenodo downloads (user provides)
│   ├── processed/
│   │   ├── isc/                       # For ISC benchmark
│   │   │   ├── <dataset_id>.rds
│   │   │   ├── <dataset_id>.rds.yaml  # Metadata + participation flags
│   │   │   └── ...
│   │   └── label_transfer/            # For label-transfer benchmark
│   │       └── <dataset_id>/
│   │           ├── query.rds
│   │           └── reference.rds
│   └── metadata/                      # Manifests, hashes for resumability
│
├── data_processing/                   # 📦 DATA PREPROCESSING
│   ├── README.md
│   ├── _targets.R                     # Targets pipeline (future)
│   ├── process_datasets.R             # Legacy script (being deprecated)
│   ├── specs_datasets.csv             # Task participation flags
│   └── config/                        # ⭐ Pipeline-specific config
│       ├── datasets.yaml              # Dataset registry + Zenodo mappings
│       └── processing_parameters.yaml # Processing parameters
│
├── ISC_benchmark/                     # 🔴 INDEPENDENT PIPELINE
│   ├── README.md                      # ISC-specific quickstart
│   ├── _targets.R                     # Targets workflow
│   ├── _targets.yaml                  # Targets configuration
│   ├── R/                             # ISC-specific functions
│   │   ├── 00_utils.R                 # Path helpers, config loading
│   │   ├── 01_isc_tasks.R             # Task definitions (6 perturbations)
│   │   ├── 02_scTypeEval_helpers.R    # Metric wrappers
│   │   └── 03_plots_tables.R          # Aggregation & plotting
│   ├── config/                        # ⭐ Pipeline-specific config
│   │   └── isc_parameters.yaml        # ISC parameters (includes global settings)
│   ├── scripts/
│   │   └── submit_hpc.sh              # SLURM submission template
│   └── results/                       # ISC-specific outputs (gitignored)
│       ├── _targets/
│       ├── task_results/
│       ├── aggregated/
│       └── figures/
│
├── label_transfer_task/               # 🟢 INDEPENDENT PIPELINE
│   ├── README.md                      # LT-specific quickstart
│   ├── _targets.R                     # Targets workflow
│   ├── _targets.yaml                  # Targets configuration
│   ├── R/                             # Label-transfer-specific functions
│   │   ├── 00_utils.R                 # Path helpers, config loading
│   │   ├── 01_classifiers.R           # Classifier definitions
│   │   ├── 02_scTypeEval_helpers.R    # Optional scTypeEval integration
│   │   └── 03_plots_tables.R          # Aggregation & plotting
│   ├── config/                        # ⭐ Pipeline-specific config
│   │   └── label_transfer_parameters.yaml # LT parameters (includes global settings)
│   ├── classifiers/
│   │   └── classifiers.R              # (Optional: separate classifier code)
│   ├── scripts/
│   │   └── submit_hpc.sh              # SLURM submission template
│   └── results/                       # LT-specific outputs (gitignored)
│       ├── _targets/
│       ├── raw_results/
│       ├── aggregated/
│       └── figures/
│
└── results/                           # ⭐ Final aggregated results (optional)
    ├── figures/
    └── tables/
```

**Key insight:** Each benchmark (🔴 ISC, 🟢 Label-Transfer) is a **stand-alone targets pipeline** that:
- Reads shared processed data from `data/processed/`
- Has its own configuration in `<benchmark>/config/`
- Produces independent results in `<benchmark>/results/`
- Can be executed in parallel (no contention)

---

## Quick Start

### 1. Install Dependencies

```bash
# Restore R environment (lock file)
Rscript -e 'renv::restore()'

# Or manually install required packages
Rscript -e '
  pkgs <- c("targets", "tarchetypes", "tidyverse", "yaml", "Seurat", 
            "SeuratObject", "Matrix")
  install.packages(pkgs)
  # remotes::install_github("carmolaab/scTypeEval")
'
```

### 2. Prepare Data

Run the data processing pipeline first (generates `data/processed/isc/` and `data/processed/label_transfer/`):

```bash
cd data_processing
Rscript -e 'targets::tar_make()'
```

See `data_processing/README.md` for details on downloading Zenodo datasets.

### 3a. Run ISC Benchmark

```bash
cd ISC_benchmark
Rscript -e 'targets::tar_make()'
```

Results → `ISC_benchmark/results/aggregated/isc_metrics_aggregated.csv`

See [ISC_benchmark/README.md](ISC_benchmark/README.md) for configuration options, HPC submission, and troubleshooting.

### 3b. Run Label-Transfer Benchmark

```bash
cd label_transfer_task
Rscript -e 'targets::tar_make()'
```

Results → `label_transfer_task/results/aggregated/label_transfer_metrics_aggregated.csv`

See [label_transfer_task/README.md](label_transfer_task/README.md) for configuration options.

### 4. (Optional) Run Both in Parallel

Since pipelines are independent, run them simultaneously:

```bash
# Terminal 1
cd ISC_benchmark && targets::tar_make()

# Terminal 2 (in parallel)
cd label_transfer_task && targets::tar_make()
```

Or submit both to HPC cluster as separate jobs.

---

## Configuration

All configuration is **YAML-based and localized to each pipeline**:

```
data_processing/config/
├── datasets.yaml                      # Dataset registry + Zenodo DOIs
└── processing_parameters.yaml         # Processing parameters (normalization, filtering, etc.)

ISC_benchmark/config/
└── isc_parameters.yaml                # All ISC settings (includes global: seed, ncores, etc.)

label_transfer_task/config/
└── label_transfer_parameters.yaml     # All LT settings (includes global: seed, ncores, etc.)
```

**Key principle:** Each pipeline is fully self-contained with all required parameters in its own `config/` directory. No shared root config required.

### Data Processing (`data_processing/config/`)

**datasets.yaml:**
```yaml
zenodo:
  doi: "10.5281/zenodo.18921437"
datasets:
  - id: Stephenson_2021_33879890
    raw_filename: "Stephenson_2021_33879890.rds"
    batch_col: "Site"
    sample_col: "Sample"
    # ...
```

**processing_parameters.yaml:**
```yaml
seed: 22
ncores: 8
processing:
  max_cells_per_sample_celltype: 200
  min_samples: 9
  max_samples: 15
  enforce_optimal_dataset: true
```

### ISC Benchmark (`ISC_benchmark/config/isc_parameters.yaml`)

```yaml
# Global settings (required for all pipelines)
seed: 22
ncores: 8

# ISC-specific settings
isc_tasks:
  run:
    missclassify: true
    Nsamples: true
    Nct: true
    NCell: true
    mergeCT: true
    SplitCelltype: true
  common:
    reduction: "umap"
    ndim: 30
    metrics:
      - silhouette
      - neighborhood_purity
    dissimilarity_methods:
      - Wasserstein
      - Pseudobulk:Euclidean

grid:
  n_replicates: 3
processing:
  min_samples: 5
  min_cells: 10
```

### Label Transfer (`label_transfer_task/config/label_transfer_parameters.yaml`)

```yaml
# Global settings (required for all pipelines)
seed: 22
ncores: 8

# Label-transfer-specific settings
classifiers:
  methods:
    - "RF"    # Random Forest
    - "SVM"   # Support Vector Machine
    - "KNN"   # k-Nearest Neighbors
  common:
    label_col: "celltype"
    pca_ndim: 30

label_transfer:
  prop_query_samples: 0.7   # Train/test split

grid:
  n_replicates: 3
```

---

## HPC Execution (SLURM)

Each pipeline includes a SLURM submission template:

### ISC Benchmark

```bash
cd ISC_benchmark
sbatch scripts/submit_hpc.sh
```

Modify `_targets.yaml` to adjust resource requests:

```yaml
controller: future
future:
  plan: future.batchtools::batchtools_slurm
  workers: 16                         # Parallel jobs
resources:
  slurm_time: "04:00:00"
  slurm_cpus_per_task: 4
  slurm_memory_gb: 32
  slurm_partition: "normal"
```

### Label-Transfer Benchmark

```bash
cd label_transfer_task
sbatch scripts/submit_hpc.sh
```

---

## Output & Results

### ISC Benchmark

- `ISC_benchmark/results/aggregated/isc_metrics_aggregated.csv` — All runs (dataset × task × rep)
- `ISC_benchmark/results/figures/` — Generated plots

### Label-Transfer Benchmark

- `label_transfer_task/results/aggregated/label_transfer_metrics_aggregated.csv` — All runs (dataset × classifier × rep)
- `label_transfer_task/results/figures/` — Generated plots

### Combining Results

Merge both aggregated tables for final paper figures:

```r
isc_results <- read.csv("ISC_benchmark/results/aggregated/isc_metrics_aggregated.csv")
lt_results <- read.csv("label_transfer_task/results/aggregated/label_transfer_metrics_aggregated.csv")
# → Merge for Figure 3 + Figure 5
```

---

## Troubleshooting

### "Could not find project root"

Ensure:
1. Working directory is within `ISC_benchmark/`, `label_transfer_task/`, or `data_processing/`
2. Project root contains `.Rproj` file and `renv.lock`

### Some datasets skipped

Check `data/processed/isc/<dataset>.rds.yaml`:

```yaml
isc: true                             # Participates in ISC benchmark
label_transfer: true                  # Participates in label-transfer
```

### Out of memory

Reduce in pipeline config:
- `replicates_per_task` / `replicates_per_classifier`
- `workers` in `_targets.yaml`
- Increase node memory in SLURM template

### Resume interrupted pipeline

Targets automatically caches completed targets:

```bash
cd ISC_benchmark
targets::tar_make()                   # Picks up where it left off
```

---

## Best Practices

### Reproducibility

- ✅ Fix random seed in pipeline-specific configs: `ISC_benchmark/config/isc_parameters.yaml` or `label_transfer_task/config/label_transfer_parameters.yaml`
- ✅ Lock R dependencies via `renv.lock`
- ✅ Log git commit hash + session info with each result
- ✅ Version YAML configs alongside results

### Development

- 📝 Edit functions in `<pipeline>/R/00_*` (write-once functions)
- 📝 Edit configs in `<pipeline>/config/` (localize parameters per pipeline)
- ✋ Never hardcode parameters in _targets.R
- ✋ Avoid nested parallelism (targets manages parallelism)

### Citation

When publishing results from this pipeline, cite:

```bibtex
@software{garnica_2026_isc_benchmark,
  title={ISC Benchmark Reproducibility},
  author={Garnica, ...},
  year={2026},
  url={https://github.com/carmolaab/ISC_benchmark_reproducibility}
}
```

---

## References

- **Targets framework:** [ropensci/targets](https://github.com/ropensci/targets)
- **scTypeEval package:** [Link to scTypeEval repo]
- **Paper:** Garnica et al., 2026

---

## License

[Specify license]

## Contact

[Specify contact info]
