# Meta_study of local and global ISC metrics

Previous results highlighted complementarity between the top task-1 metric and the top task-2 metric:

- task-1-oriented metric: 2-label silhouette on pseudobulk cosine dissimilarity
- task-2-oriented metric: silhouette on RCM dissimilarity

In the manuscript, these are referred to as:

- ISCG (global ISC): cosine-based ISC, which captures global cell type compactness (cross-sample similarity within a cell type relative to all other cell types)
- ISCL (local ISC): RCM-based ISC, which captures cell-type inter-sample replicability (how confidently a cell type profile maps across samples, especially for closely related and hard-to-distinguish cell types)

To characterize real-world ISC ranges and derive empirical score-interpretation guidelines, the Meta-study analysis evaluates cell type consistency across nine scRNA-seq datasets spanning diverse tissues and study designs.

## Outputs

The pipeline writes scTypeEval objects to:

- output/scTypeEval_objs

Each output file follows this naming scheme:

- <dataset_ds_name>__<annotation_column>.rds

These files are consumed by Figures_notebooks/5_Figure.Rmd.

## Configuration

- config/meta_study_parameters.yaml: paths, execution settings, dissimilarity methods.
- config/datasets_metadata.yaml: dataset metadata and annotation columns.


