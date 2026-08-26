#!/usr/bin/env python3
"""Run CellHint harmonization and export per-cell-type scores.

This script is intended to mirror the per-cluster output style used for SCCAF,
but with CellHint's harmonization signal across datasets/samples.

Output score interpretation:
- dominant_group_fraction near 1.0: the cell type maps mostly to one
  harmonized group (high consistency).
- lower values: the cell type is split across multiple harmonized groups.
"""

import argparse
import importlib
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

# This script name (cellhint.py) shadows the package name.
# Remove this directory from sys.path before importing the external package.
_SCRIPT_DIR = str(Path(__file__).resolve().parent)
if _SCRIPT_DIR in sys.path:
	sys.path.remove(_SCRIPT_DIR)
cellhint_pkg = importlib.import_module("cellhint")


def _pick_column(adata, requested_key, preferred, label):
	"""Pick an obs column by explicit key or preferred fallback list."""
	if requested_key:
		if requested_key not in adata.obs.columns:
			raise ValueError(
				f"{label} column '{requested_key}' not found in adata.obs. "
				f"Available columns: {list(adata.obs.columns)}"
			)
		return requested_key

	for key in preferred:
		if key in adata.obs.columns:
			return key

	raise ValueError(
		f"No {label} column provided and none of the default columns were found. "
		f"Available columns: {list(adata.obs.columns)}"
	)


def _ensure_pca(adata, n_top_genes=3000, n_pcs=50):
	"""Ensure `X_pca` exists in obsm; compute it from normalized data if needed."""
	if "X_pca" in adata.obsm:
		return adata

	ad = adata.copy()
	sc.pp.normalize_total(ad, target_sum=1e4)
	sc.pp.log1p(ad)

	n_top = min(int(n_top_genes), ad.n_vars)
	sc.pp.highly_variable_genes(ad, n_top_genes=n_top, flavor="cell_ranger")
	if "highly_variable" in ad.var.columns and int(ad.var["highly_variable"].sum()) > 10:
		ad = ad[:, ad.var["highly_variable"]].copy()

	max_pcs = max(2, min(int(n_pcs), ad.n_obs - 1, ad.n_vars - 1))
	sc.pp.pca(ad, n_comps=max_pcs)
	return ad


def run_cellhint_per_cluster(
	h5ad_path,
	output_csv,
	cluster_key=None,
	dataset_key=None,
	random_state=2,
):
	"""Run CellHint harmonization and export per-cell-type scores."""
	adata = sc.read_h5ad(h5ad_path)

	cluster_col = _pick_column(
		adata,
		cluster_key,
		preferred=[
			"OriginalAnnotationLevel2",
			"OriginalAnnotationLevel1",
			"cell_type",
			"celltype",
			"annotation",
			"louvain",
			"leiden",
		],
		label="cluster",
	)
	dataset_col = _pick_column(
		adata,
		dataset_key,
		preferred=["sample", "dataset", "batch", "donor", "patient"],
		label="dataset",
	)

	adata = _ensure_pca(adata)

	# Harmonize cell types across datasets/samples using PCA representation.
	alignment = cellhint_pkg.harmonize(
		adata,
		dataset=dataset_col,
		cell_type=cluster_col,
		use_rep="X_pca",
		random_state=int(random_state),
	)

	# reannotation has per-cell assignment to harmonized groups.
	reann = alignment.reannotation[["dataset", "cell_type", "group"]].copy()
	reann["dataset"] = reann["dataset"].astype(str)
	reann["cell_type"] = reann["cell_type"].astype(str)
	reann["group"] = reann["group"].astype(str)

	rows = []
	input_path = str(Path(h5ad_path).resolve())
	n_total_datasets = int(reann["dataset"].nunique())

	for cell_type, block in reann.groupby("cell_type"):
		n_cells = int(len(block))
		group_counts = block["group"].value_counts()
		dominant_group = str(group_counts.index[0])
		dominant_count = int(group_counts.iloc[0])
		dominant_fraction = float(dominant_count / n_cells)

		# Entropy gives how spread labels are across harmonized groups.
		probs = (group_counts / n_cells).to_numpy(dtype=float)
		entropy = float(max(0.0, -(probs * np.log2(probs + 1e-12)).sum()))
		max_entropy = float(np.log2(max(len(group_counts), 2)))
		normalized_entropy = float(np.clip(entropy / max_entropy, 0.0, 1.0))

		by_dataset = (
			block.groupby("dataset")["group"].agg(lambda x: x.value_counts().index[0])
		)
		datasets_with_dominant = int((by_dataset == dominant_group).sum())

		rows.append(
			{
				"input_h5ad": input_path,
				"dataset_key": dataset_col,
				"cluster_key": cluster_col,
				"celltype": cell_type,
				"score": dominant_fraction,
				"cell_type": cell_type,
				"n_cells": n_cells,
				"n_harmonized_groups": int(len(group_counts)),
				"dominant_group": dominant_group,
				"dominant_group_fraction": dominant_fraction,
				"datasets_total": n_total_datasets,
				"datasets_with_dominant_group": datasets_with_dominant,
				"dominant_group_dataset_fraction": float(datasets_with_dominant / n_total_datasets),
				"group_entropy": entropy,
				"group_entropy_normalized": normalized_entropy,
				"consistency_score": dominant_fraction,
			}
		)

	out_df = pd.DataFrame(rows).sort_values(
		["consistency_score", "n_cells"], ascending=[False, False]
	)

	output_path = Path(output_csv)
	output_path.parent.mkdir(parents=True, exist_ok=True)
	out_df.to_csv(output_path, index=False)
	return output_path


def main():
	parser = argparse.ArgumentParser(
		description="Run CellHint and export per-cell-type consistency-like scores."
	)
	parser.add_argument("--input", required=True, help="Path to input .h5ad")
	parser.add_argument(
		"--output",
		default=None,
		help="Output CSV path (default: <input_stem>_cellhint_per_cluster.csv)",
	)
	parser.add_argument(
		"--cluster-key",
		default=None,
		help="obs column with cell type / cluster labels (auto-detect if omitted)",
	)
	parser.add_argument(
		"--dataset-key",
		default=None,
		help="obs column defining samples/datasets/batches (auto-detect if omitted)",
	)
	parser.add_argument(
		"--random-state",
		type=int,
		default=2,
		help="Random seed passed to cellhint.harmonize (default: 2)",
	)
	args = parser.parse_args()

	input_path = Path(args.input)
	output_csv = args.output or str(
		input_path.with_name(f"{input_path.stem}_cellhint_per_cluster.csv")
	)

	output_path = run_cellhint_per_cluster(
		h5ad_path=str(input_path),
		output_csv=output_csv,
		cluster_key=args.cluster_key,
		dataset_key=args.dataset_key,
		random_state=args.random_state,
	)
	print(f"Saved CellHint per-cell-type scores to: {output_path}")


if __name__ == "__main__":
	main()