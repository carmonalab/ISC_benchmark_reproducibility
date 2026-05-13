# Meta_study Targets Pipeline

This pipeline reproduces the scTypeEval object generation logic from MetaStudy_consistency.Rmd using an incremental targets workflow.

## Outputs

The pipeline writes scTypeEval objects to:

- output/scTypeEval_objs

Each output file follows this naming scheme:

- <dataset_ds_name>__<annotation_column>.rds

These files are consumed by Figures_notebooks/5_Figure.Rmd.

## Configuration

- config/meta_study_parameters.yaml: paths, execution settings, dissimilarity methods.
- config/datasets_metadata.yaml: dataset metadata and annotation columns.

## Run

From Meta_study:

- targets::tar_make()

Or with shell script:

- ./scripts/run_targets.sh

## Incremental behavior

- One target is created per dataset x annotation column.
- Existing output files are skipped when execution.skip_existing = true.
- Changing YAML config invalidates only affected downstream targets.

## Notes

- If package scTypeEval is unavailable, the pipeline falls back to sourcing R scripts from paths.scTypeEval_source_dir.
- The default blacklist is loaded from paths.default_black_list_rdata.
