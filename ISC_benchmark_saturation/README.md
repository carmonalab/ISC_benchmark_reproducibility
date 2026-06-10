# ISC Benchmark Saturation Analysis Pipeline

Ranking stability analysis that evaluates method robustness across progressively larger subsets of datasets.

## Scope

This pipeline assesses the stability of method rankings (consistency metrics) as a function of dataset selection:
- Evaluates ranking stability by computing Spearman correlations between method rankings from disjoint dataset subsets
- Measures ranking convergence at increasing dataset subset sizes
- Fits trend curves (generalized additive models) to predict ranking stability at larger dataset numbers
- Identifies whether top-performing methods remain stable across different dataset compositions

Main entrypoints:
- `_targets.R`
- `master_job.sh`
- `scripts/submit_saturation_hpc.sh`

## Inputs

Required:
- ISC benchmark results from `../results/isc_benchmark/all_results_merged.rds`
- Pipeline parameters in `config/saturation_parameters.yaml`
