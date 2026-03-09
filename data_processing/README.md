# Dataset processing

This directory contains the **reproducible dataset preprocessing** used for the manuscript.

## Inputs

- `data_processing/specs_datasets.csv` — dataset/condition/batch specifications and which benchmarks apply
- `config/datasets.yaml` — maps each dataset reference to a Zenodo file and metadata column mappings
- `config/benchmark_parameters.yaml` — downsampling + split parameters
- `zenodo/raw/` — raw Zenodo `.rds` files (downloaded separately)

## Zenodo datasets

This repository does **not** ship the raw single-cell datasets.

- Zenodo DOI: `10.5281/zenodo.18921437`

### Download

Option A (manual):
1. Download all `.rds` files from the Zenodo record.
2. Move them into `zenodo/raw/`.

Option B (scripted):
- Populate direct file URLs in `config/datasets.yaml`, then run:
	- `Rscript data_processing/process_datasets.R --download true`

The processing script will skip files that already exist locally.

## Output

- `data/processed/isc/*.rds` — processed Seurat objects (one per dataset-batch-condition)
- `data/processed/label_transfer/<dataset_id>/{query.rds,reference.rds}` — query/reference splits for label transfer

## Run

From the repository root:

- `Rscript data_processing/process_datasets.R`
