#!/usr/bin/env python3
"""Run popV and export per-cell-type likelihood from classifier agreement.

This script uses popV's classifier-agreement score (number of classifiers
agreeing on the final prediction) as the certainty signal. It then aggregates
that signal per query cell type into a 0-1 likelihood.

Default behavior uses the pretrained Hub model `popV/tabula_sapiens_All_Cells`
in `fast` mode, matching popV's fast annotation workflow.
"""

import argparse
import importlib
import json
import sys
import warnings
from pathlib import Path

import anndata as ad
import pandas as pd

# This script name shadows the package name.
_SCRIPT_DIR = str(Path(__file__).resolve().parent)
if _SCRIPT_DIR in sys.path:
	sys.path.remove(_SCRIPT_DIR)

popv_pkg = importlib.import_module("popv")
popv_preprocessing = importlib.import_module("popv.preprocessing")
popv_annotation = importlib.import_module("popv.annotation")
popv_hub = importlib.import_module("popv.hub")


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
	fast_methods = set(popv_annotation.algorithms_nt.FAST_ALGORITHMS)
	valid = set(popv_annotation.algorithms_nt.ALL_ALGORITHMS)

	if methods_arg is None or methods_arg.strip() == "":
		if prediction_mode == "fast":
			return sorted(fast_methods)
		return list(popv_annotation.algorithms_nt.CURRENT_ALGORITHMS)

	methods = [m.strip() for m in methods_arg.split(",") if m.strip()]
	invalid = [m for m in methods if m not in valid]
	if invalid:
		raise ValueError(
			f"Invalid method(s): {invalid}. Valid methods include: {sorted(valid)}"
		)

	if prediction_mode == "fast":
		non_fast = [m for m in methods if m not in fast_methods]
		if non_fast:
			raise ValueError(
				"In fast mode, only pretrained-compatible methods are allowed. "
				f"Unsupported in fast mode: {non_fast}. "
				f"Use only: {sorted(fast_methods)}"
			)
	return methods


def _to_numeric(series):
	return pd.to_numeric(series.astype(str), errors="coerce")


def _get_query_obs(adata):
	if "_dataset" in adata.obs.columns:
		query_mask = adata.obs["_dataset"].astype(str).eq("query")
		if query_mask.any():
			return adata.obs.loc[query_mask].copy()
	return adata.obs.copy()


def _build_summary(
	query_obs,
	query_h5ad,
	output_csv,
	prediction_mode,
	selected_methods,
	query_celltype_key=None,
	model_source="",
	reference_h5ad="",
	ref_labels_key="",
	ref_batch_key="",
	query_batch_key="",
):
	if "popv_prediction_score" not in query_obs.columns:
		raise ValueError(
			"popv_prediction_score not found in popV outputs. "
			"Check popV run configuration."
		)

	n_methods = len(selected_methods)
	if n_methods <= 0:
		raise ValueError("No prediction methods available after selection.")

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
	summary["celltype"] = summary["query_cell_type"].astype(str)
	summary["score"] = summary["agreement_likelihood_mean"].astype(float)

	summary["n_methods"] = int(n_methods)
	summary["query_h5ad"] = str(Path(query_h5ad).resolve())
	summary["reference_h5ad"] = reference_h5ad
	summary["ref_labels_key"] = ref_labels_key
	summary["ref_batch_key"] = ref_batch_key
	summary["query_batch_key"] = query_batch_key
	summary["query_celltype_key"] = query_celltype_key if query_celltype_key is not None else "popv_prediction"
	summary["prediction_mode"] = prediction_mode
	summary["methods"] = ",".join(selected_methods)
	summary["model_source"] = model_source

	summary = summary.sort_values(
		["agreement_likelihood_mean", "n_cells"],
		ascending=[False, False],
	)

	output_path = Path(output_csv)
	output_path.parent.mkdir(parents=True, exist_ok=True)
	summary.to_csv(output_path, index=False)
	return output_path


def _load_reference_genes_from_hub_model(hub_model):
	"""Load reference gene IDs expected by the pretrained Hub model."""
	preprocessing_json = Path(hub_model.local_dir) / "preprocessing.json"
	if not preprocessing_json.exists():
		raise ValueError(
			f"Missing preprocessing config in Hub model cache: {preprocessing_json}."
		)

	with open(preprocessing_json) as handle:
		data = json.load(handle)

	genes = data.get("gene_names")
	if not genes:
		raise ValueError(
			f"No gene_names found in Hub preprocessing config: {preprocessing_json}."
		)
	return set(genes)


def _overlap_ratio(query_var_names, reference_genes):
	if not reference_genes:
		return 0.0
	query_genes = set(map(str, query_var_names))
	return len(query_genes.intersection(reference_genes)) / float(len(reference_genes))


def _annotate_hub_with_gene_mapping(
	hub_model,
	query_adata,
	query_batch_key,
	prediction_mode,
	selected_methods,
):
	"""Annotate query data with automated gene-symbol mapping strategy.

	Strategy:
	1) Map using feature_name and compute overlap to pretrained reference genes.
	2) If overlap is <50%, retry mapping using gene_symbol.
	"""
	annotate_kwargs = {
		"query_batch_key": query_batch_key,
		"save_path": "tmp/popv_output",
		"prediction_mode": prediction_mode,
		"methods": selected_methods,
	}
	reference_genes = _load_reference_genes_from_hub_model(hub_model)

	# Try feature_name mapping first.
	feature_mapped = hub_model.map_genes(
		query_adata.copy(),
		gene_symbols="feature_name",
		organism=hub_model.metadata.organism,
	)
	feature_overlap = _overlap_ratio(feature_mapped.var_names, reference_genes)

	if feature_overlap >= 0.5:
		return hub_model.annotate_data(
			query_adata=feature_mapped,
			**annotate_kwargs,
		)

	warnings.warn(
		"Gene overlap after feature_name mapping is below 50% "
		f"({feature_overlap:.2%}); retrying with gene_symbol mapping.",
		UserWarning,
		stacklevel=2,
	)

	gene_symbol_mapped = hub_model.map_genes(
		query_adata.copy(),
		gene_symbols="gene_symbol",
		organism=hub_model.metadata.organism,
	)
	gene_symbol_overlap = _overlap_ratio(gene_symbol_mapped.var_names, reference_genes)

	if gene_symbol_overlap < 0.5:
		raise ValueError(
			"Low overlap between mapped query genes and pretrained reference genes. "
			f"feature_name overlap={feature_overlap:.2%}, "
			f"gene_symbol overlap={gene_symbol_overlap:.2%}. "
			"Please verify query gene identifiers or switch to query/reference mode."
		)

	return hub_model.annotate_data(
		query_adata=gene_symbol_mapped,
		**annotate_kwargs,
	)


def run_popv_pretrained_hub(
	query_h5ad,
	output_csv,
	pretrained_model_repo,
	query_batch_key=None,
	query_celltype_key=None,
	prediction_mode="fast",
	methods=None,
	hub_cache_dir=None,
):
	"""Run popV using a pretrained Hub model and export query certainty summary."""
	query_adata = ad.read_h5ad(query_h5ad)

	if query_batch_key is not None:
		_require_obs_column(query_adata, query_batch_key, "query")

	hub_model = popv_hub.HubModel.pull_from_huggingface_hub(
		pretrained_model_repo,
		cache_dir=hub_cache_dir,
	)

	selected_methods = _resolve_methods(methods, prediction_mode)
	annotated_adata = _annotate_hub_with_gene_mapping(
		hub_model=hub_model,
		query_adata=query_adata,
		query_batch_key=query_batch_key,
		prediction_mode=prediction_mode,
		selected_methods=selected_methods,
	)

	query_obs = _get_query_obs(annotated_adata)
	return _build_summary(
		query_obs=query_obs,
		query_h5ad=query_h5ad,
		output_csv=output_csv,
		prediction_mode=prediction_mode,
		selected_methods=selected_methods,
		query_celltype_key=query_celltype_key,
		model_source=pretrained_model_repo,
		query_batch_key=query_batch_key if query_batch_key is not None else "",
	)


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
	ref_batch_key = _ensure_batch_column(ref_adata, ref_batch_key, "reference", "_popv_ref_batch")
	query_batch_key = _ensure_batch_column(query_adata, query_batch_key, "query", "_popv_query_batch")

	selected_methods = _resolve_methods(methods, prediction_mode)
	if not selected_methods:
		raise ValueError(
			"No usable popV methods remain after method filtering. "
			"Provide explicit methods with --methods."
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

	query_obs = _get_query_obs(adata)
	return _build_summary(
		query_obs=query_obs,
		query_h5ad=query_h5ad,
		output_csv=output_csv,
		prediction_mode=prediction_mode,
		selected_methods=selected_methods,
		query_celltype_key=query_celltype_key,
		model_source="query_reference",
		reference_h5ad=str(Path(reference_h5ad).resolve()),
		ref_labels_key=ref_labels_key,
		ref_batch_key=ref_batch_key,
		query_batch_key=query_batch_key if query_batch_key is not None else "",
	)


def main():
	parser = argparse.ArgumentParser(
		description=(
			"Run popV and export per-query-cell-type likelihood scores. "
			"Defaults to Tabula Sapiens pretrained Hub model in fast mode."
		)
	)
	parser.add_argument("--query", required=True, help="Path to query .h5ad")
	parser.add_argument(
		"--reference",
		default=None,
		help="Optional path to reference .h5ad. If omitted, Hub pretrained model is used.",
	)
	parser.add_argument(
		"--output",
		default=None,
		help="Output CSV path (default: <query_stem>_popv_per_celltype_scores.csv)",
	)
	parser.add_argument(
		"--ref-labels-key",
		default=None,
		help="Reference obs column with cell-type labels (required when --reference is provided)",
	)
	parser.add_argument(
		"--ref-batch-key",
		default=None,
		help="Optional reference obs column with batch/sample annotation (query/reference mode only).",
	)
	parser.add_argument(
		"--query-batch-key",
		default=None,
		help="Optional query obs column with batch/sample annotation.",
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
		help="popV prediction mode (default: fast).",
	)
	parser.add_argument(
		"--methods",
		default=None,
		help=(
			"Comma-separated popV methods to use. "
			"In fast mode, only FAST methods are allowed."
		),
	)
	parser.add_argument(
		"--pretrained-model-repo",
		default="popV/tabula_sapiens_All_Cells",
		help="HuggingFace Hub model repo for pretrained popV model (default: popV/tabula_sapiens_All_Cells).",
	)
	parser.add_argument(
		"--hub-cache-dir",
		default=None,
		help="Optional cache directory for downloaded Hub model artifacts.",
	)
	parser.add_argument(
		"--n-samples-per-label",
		type=int,
		default=300,
		help="Reference subsampling per label during preprocessing (query/reference mode only).",
	)
	parser.add_argument(
		"--hvg",
		type=int,
		default=4000,
		help="Number of HVGs for preprocessing (query/reference mode only).",
	)
	parser.add_argument(
		"--unknown-celltype-label",
		default="unknown",
		help="Unknown label used by popV preprocessing (query/reference mode only).",
	)
	parser.add_argument(
		"--save-path-trained-models",
		default="tmp/popv_models",
		help="Directory for popV temporary/pretrained model artifacts (query/reference mode only).",
	)

	args = parser.parse_args()
	query_path = Path(args.query)
	output_csv = args.output or str(
		query_path.with_name(f"{query_path.stem}_popv_per_celltype_scores.csv")
	)

	if args.reference is None:
		if args.prediction_mode != "fast":
			raise ValueError(
				"Pretrained Hub workflow without --reference supports only --prediction-mode fast."
			)
		out_path = run_popv_pretrained_hub(
			query_h5ad=args.query,
			output_csv=output_csv,
			pretrained_model_repo=args.pretrained_model_repo,
			query_batch_key=args.query_batch_key,
			query_celltype_key=args.query_celltype_key,
			prediction_mode=args.prediction_mode,
			methods=args.methods,
			hub_cache_dir=args.hub_cache_dir,
		)
	else:
		if args.ref_labels_key is None:
			raise ValueError("--ref-labels-key is required when --reference is provided.")
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
