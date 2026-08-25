#!/usr/bin/env python3
"""Run popV with query/reference h5ad files and export per-cell-type likelihood.

This script uses popV's classifier-agreement score (number of classifiers
agreeing on the final prediction) as the certainty signal. It then aggregates
that signal per query cell type into a 0-1 likelihood.
"""

import argparse
import importlib
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

# This script name shadows the package name.
_SCRIPT_DIR = str(Path(__file__).resolve().parent)
if _SCRIPT_DIR in sys.path:
	sys.path.remove(_SCRIPT_DIR)

popv_pkg = importlib.import_module("popv")
popv_preprocessing = importlib.import_module("popv.preprocessing")
popv_annotation = importlib.import_module("popv.annotation")


def _require_obs_column(adata, key, dataset_name):
	if key not in adata.obs.columns:
		raise ValueError(
			f"{dataset_name} is missing required obs column '{key}'. "
			f"Available columns: {list(adata.obs.columns)}"
		)


def _ensure_batch_column(adata, batch_key, dataset_name, default_name):
	"""Return an obs column name for batch annotation, creating a dummy one if needed."""
	if batch_key is not None:
		_require_obs_column(adata, batch_key, dataset_name)
		return batch_key

	adata.obs[default_name] = "batch_1"
	return default_name


def _resolve_methods(methods_arg, prediction_mode):
	if methods_arg is None or methods_arg.strip() == "":
		if prediction_mode == "fast":
			return list(popv_annotation.algorithms_nt.FAST_ALGORITHMS)
		return list(popv_annotation.algorithms_nt.CURRENT_ALGORITHMS)

	methods = [m.strip() for m in methods_arg.split(",") if m.strip()]
	valid = set(popv_annotation.algorithms_nt.ALL_ALGORITHMS)
	invalid = [m for m in methods if m not in valid]
	if invalid:
		raise ValueError(
			f"Invalid method(s): {invalid}. Valid methods include: {sorted(valid)}"
		)
	return methods


def _filter_methods_for_batch_support(methods, has_batch_support):
	if has_batch_support:
		return methods

	batch_methods = {"KNN_BBKNN", "KNN_HARMONY"}
	filtered = [method for method in methods if method not in batch_methods]
	removed = sorted(set(methods) - set(filtered))
	if removed:
		print(
			"Skipping batch-aware popV methods because no batch keys were provided: "
			+ ", ".join(removed)
		)
	return filtered


def _to_numeric(series):
	return pd.to_numeric(series.astype(str), errors="coerce")


def run_popv_query_reference(
	query_h5ad,
	reference_h5ad,
	output_csv,
	ref_labels_key,
	ref_batch_key,
	query_batch_key=None,
	query_celltype_key=None,
	prediction_mode="fast",
	methods=None,
	n_samples_per_label=300,
	hvg=4000,
	unknown_celltype_label="unknown",
	save_path_trained_models="tmp/popv_models",
):
	"""Run popV and export per-query-cell-type certainty summary."""
	query_adata = ad.read_h5ad(query_h5ad)
	ref_adata = ad.read_h5ad(reference_h5ad)

	_require_obs_column(ref_adata, ref_labels_key, "reference")
	has_batch_support = ref_batch_key is not None or query_batch_key is not None
	ref_batch_key = _ensure_batch_column(ref_adata, ref_batch_key, "reference", "_popv_ref_batch")
	query_batch_key = _ensure_batch_column(query_adata, query_batch_key, "query", "_popv_query_batch")

	selected_methods = _resolve_methods(methods, prediction_mode)
	selected_methods = _filter_methods_for_batch_support(selected_methods, has_batch_support)
	if not selected_methods:
		raise ValueError(
			"No usable popV methods remain after removing batch-aware methods. "
			"Provide batch annotations or choose explicit non-batch methods."
		)

	processor = popv_preprocessing.Process_Query(
		query_adata=query_adata,
		ref_adata=ref_adata,
		ref_labels_key=ref_labels_key,
		ref_batch_key=ref_batch_key,
		query_batch_key=query_batch_key,
		# Disable ontology dependency for generic reproducibility benchmark usage.
		cl_obo_folder=False,
		prediction_mode=prediction_mode,
		unknown_celltype_label=unknown_celltype_label,
		n_samples_per_label=n_samples_per_label,
		save_path_trained_models=save_path_trained_models,
		hvg=hvg,
	)
	adata = processor.adata

	popv_annotation.annotate_data(adata, methods=selected_methods, save_path=None)

	query_mask = adata.obs["_dataset"].astype(str).eq("query")
	query_obs = adata.obs.loc[query_mask].copy()

	# In ontology-disabled mode, popv_prediction_score equals majority-vote agreement.
	if "popv_prediction_score" not in query_obs.columns:
		raise ValueError(
			"popv_prediction_score not found in popV outputs. "
			"Check popV run configuration."
		)

	n_methods = len(adata.uns.get("prediction_keys_seen", []))
	if n_methods <= 0:
		raise ValueError("No prediction methods found in popV output (prediction_keys_seen is empty).")

	# Agreement count per cell, then normalized likelihood in [0,1].
	query_obs["agreement_count"] = _to_numeric(query_obs["popv_prediction_score"])
	query_obs["agreement_likelihood_0_1"] = query_obs["agreement_count"] / float(n_methods)

	if query_celltype_key is not None:
		_require_obs_column(query_obs, query_celltype_key, "query results")
		group_col = query_celltype_key
	else:
		# Default to popV prediction as query annotation if user does not provide one.
		group_col = "popv_prediction"

	summary = (
		query_obs.groupby(group_col, dropna=False)
		.agg(
			n_cells=(group_col, "size"),
			agreement_likelihood_mean=("agreement_likelihood_0_1", "mean"),
			agreement_likelihood_median=("agreement_likelihood_0_1", "median"),
			agreement_likelihood_sd=("agreement_likelihood_0_1", "std"),
			agreement_count_mean=("agreement_count", "mean"),
		)
		.reset_index()
		.rename(columns={group_col: "query_cell_type"})
	)

	summary["n_methods"] = int(n_methods)
	summary["query_h5ad"] = str(Path(query_h5ad).resolve())
	summary["reference_h5ad"] = str(Path(reference_h5ad).resolve())
	summary["ref_labels_key"] = ref_labels_key
	summary["ref_batch_key"] = ref_batch_key
	summary["query_batch_key"] = query_batch_key if query_batch_key is not None else ""
	summary["query_celltype_key"] = query_celltype_key if query_celltype_key is not None else "popv_prediction"
	summary["prediction_mode"] = prediction_mode
	summary["methods"] = ",".join(selected_methods)

	summary = summary.sort_values(
		["agreement_likelihood_mean", "n_cells"],
		ascending=[False, False],
	)

	output_path = Path(output_csv)
	output_path.parent.mkdir(parents=True, exist_ok=True)
	summary.to_csv(output_path, index=False)
	return output_path


def main():
	parser = argparse.ArgumentParser(
		description="Run popV on query/reference h5ad and export per-query-cell-type likelihood scores."
	)
	parser.add_argument("--query", required=True, help="Path to query .h5ad")
	parser.add_argument("--reference", required=True, help="Path to reference .h5ad")
	parser.add_argument(
		"--output",
		default=None,
		help="Output CSV path (default: <query_stem>_popv_per_celltype_scores.csv)",
	)
	parser.add_argument(
		"--ref-labels-key",
		required=True,
		help="Reference obs column with cell-type labels",
	)
	parser.add_argument(
		"--ref-batch-key",
		default=None,
		help="Optional reference obs column with batch/sample annotation. If omitted, a dummy batch is created.",
	)
	parser.add_argument(
		"--query-batch-key",
		default=None,
		help="Optional query obs column with batch/sample annotation",
	)
	parser.add_argument(
		"--query-celltype-key",
		default=None,
		help="Optional query obs column to aggregate scores by. Default: popv_prediction",
	)
	parser.add_argument(
		"--prediction-mode",
		choices=["fast", "retrain", "inference"],
		default="fast",
		help="popV prediction mode (default: fast)",
	)
	parser.add_argument(
		"--methods",
		default=None,
		help="Comma-separated popV methods to use. Default: FAST methods in fast mode, CURRENT otherwise.",
	)
	parser.add_argument(
		"--n-samples-per-label",
		type=int,
		default=300,
		help="Reference subsampling per label during preprocessing (default: 300)",
	)
	parser.add_argument(
		"--hvg",
		type=int,
		default=4000,
		help="Number of HVGs for preprocessing (default: 4000)",
	)
	parser.add_argument(
		"--unknown-celltype-label",
		default="unknown",
		help="Unknown label used by popV preprocessing (default: unknown)",
	)
	parser.add_argument(
		"--save-path-trained-models",
		default="tmp/popv_models",
		help="Directory for popV temporary/pretrained model artifacts",
	)

	args = parser.parse_args()
	query_path = Path(args.query)
	output_csv = args.output or str(
		query_path.with_name(f"{query_path.stem}_popv_per_celltype_scores.csv")
	)

	out_path = run_popv_query_reference(
		query_h5ad=args.query,
		reference_h5ad=args.reference,
		output_csv=output_csv,
		ref_labels_key=args.ref_labels_key,
		ref_batch_key=args.ref_batch_key,
		query_batch_key=args.query_batch_key,
		query_celltype_key=args.query_celltype_key,
		prediction_mode=args.prediction_mode,
		methods=args.methods,
		n_samples_per_label=args.n_samples_per_label,
		hvg=args.hvg,
		unknown_celltype_label=args.unknown_celltype_label,
		save_path_trained_models=args.save_path_trained_models,
	)
	print(f"Saved popV per-query-cell-type likelihood scores to: {out_path}")


if __name__ == "__main__":
	main()
