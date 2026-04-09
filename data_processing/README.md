# Dataset Processing Pipeline

Reproduces the core data preprocessing from `Consistency_metrics_benchmark/datasets/proc_data.Rmd`.

**Processes 6 datasets**: downsamples cells, splits by batch × condition, applies quality filters.

---

## Quick Start

### 1. Download Raw Data

Download 6 raw Seurat objects from Zenodo (10.5281/zenodo.18921437) to `data/raw/`:
```
JoaI_2022_35773407_Nofilt_whole.rds
StephensonE_2021_33879890_preprocessed.rds
BCC_all_LY_clean_annotated.rds
BCC_CG_all_annotated.rds
LUNG_b351804c-293e-4aeb-9c4c-043db67f4540.rds
ICB_8d918bdd-ab11-4c83-9de0-93640aeb8e20.rds
```

### 2. Run Processing

From repository root:
```bash
Rscript data_processing/process_datasets.R
```

### 3. Verify Output

```bash
ls data/processed/*.rds | wc -l   # Should be 31 files
```

---

## Configuration

**`config/processing_parameters.yaml`**:
```yaml
seed: 22              # Reproducibility seed (matches proc_data.Rmd)
n_cores: 8            # Parallel cores

downsampling:
  max_cells_per_sample_celltype: 200

sample_filtering:
  min_samples: 9      # Minimum samples per split
  max_samples: 15     # Maximum samples per split
  min_cells_per_sample: 200
```

**`config/core_datasets.yaml`**: Registry of 6 core datasets with metadata mappings

---

## Pipeline Steps

For each dataset:

1. **Load** raw Seurat object
2. **Exclude** low-quality cell types (dataset-specific)
3. **Downsample** cells (max 200 per cell type per sample, seed=22)
4. **Split** by batch + condition
5. **Filter** each split with `optimal_dataset()` (validates sample/cell requirements)
6. **Save** as `PREFIX_batch_condition.rds`

**Output**: 31 files matching reference implementation

---

## Core Datasets Overview

| Dataset | Reference | Splits | Notes |
|---------|-----------|--------|-------|
| JoaI | Joanito et al. 2022 | 7 | Colorectal cancer (dataset × tumor/normal) |
| Mitchel | Stephenson et al. 2021 | 5 | COVID PBMC (site × status) |
| BCC | Ganier & Yerly | 2 | Basal cell carcinoma (dataset only) |
| LungAtlas | Sikkema et al. 2023 | 8 | Healthy lung (dataset × tissue) |
| ICBAtlas | Gondal et al. 2025 | 9 | Immunotherapy cohort (study × pre/post) |

---


## Key Functions

**`process_dataset()`**: Main processing function
- Loads, filters, downsamples, splits, and saves one dataset
- Called iteratively for each dataset in registry

**`optimal_dataset()`**: Quality control filter
- Validates sample count (min/max bounds)
- Selects top samples by cell type diversity + cell count
- Returns NULL if split fails QC

---

