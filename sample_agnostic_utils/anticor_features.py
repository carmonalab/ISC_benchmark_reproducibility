#!/usr/bin/env python3
"""Run anti-correlation feature discovery per cell type and export CSV.

For each cell type/cluster in an h5ad object, this script runs
`get_anti_cor_genes` on the corresponding subset of cells and reports:
- number of selected anti-correlated genes
- an asymptotic 0-1 score based on anti-correlated genes relative to
	expressed genes in that cell type
"""

import argparse
import importlib
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

# This script name shadows the package name. Remove this directory from sys.path
# so we import the installed package, not this file.
_SCRIPT_DIR = str(Path(__file__).resolve().parent)
if _SCRIPT_DIR in sys.path:
	sys.path.remove(_SCRIPT_DIR)
anticor_pkg = importlib.import_module("anticor_features")


def _pick_cluster_key(adata, requested_key=None):
	"""Pick a cluster column from adata.obs."""
	if requested_key:
		if requested_key not in adata.obs.columns:
			raise ValueError(
				f"Cluster column '{requested_key}' not found in adata.obs. "
				f"Available columns: {list(adata.obs.columns)}"
			)
		return requested_key

	preferred = [
		"OriginalAnnotationLevel2",
		"OriginalAnnotationLevel1",
		"cell_type",
		"celltype",
		"annotation",
		"louvain",
		"leiden",
		"seurat_clusters",
	]
	for key in preferred:
		if key in adata.obs.columns:
			return key

	raise ValueError(
		"No cluster column provided and no default cluster columns were found. "
		f"Available columns: {list(adata.obs.columns)}"
	)


def _to_gene_by_cell_matrix(x):
	"""Convert AnnData X (cells x genes) to genes x cells format expected by package."""
	if sparse.issparse(x):
		return x.T.tocsr()
	return np.asarray(x).T


def _scale_minmax(values):
	"""Deprecated helper retained for backward compatibility."""
	arr = np.asarray(values, dtype=float)
	return arr


def _count_expressed_genes(x_sub):
	"""Count genes with expression > 0 in at least one cell in subset."""
	if sparse.issparse(x_sub):
		return int(np.sum(np.asarray(x_sub.getnnz(axis=0) > 0)))
	x_arr = np.asarray(x_sub)
	return int(np.sum(np.any(x_arr > 0, axis=0)))


def _expressed_gene_mask(x_sub):
	"""Boolean mask for genes expressed (>0) in at least one cell of subset."""
	if sparse.issparse(x_sub):
		return np.asarray(x_sub.getnnz(axis=0) > 0).ravel()
	x_arr = np.asarray(x_sub)
	return np.any(x_arr > 0, axis=0)


def _asymptotic_score(anticor_gene_count, n_features_expressed, score_k):
	"""Asymptotic score in [0,1], with score=1 when anticor_gene_count=0.

	Let r = anticor_gene_count / n_features_expressed.
	Score is exp(-k * r), so it drops rapidly near r=0 and flattens later.
	"""
	if n_features_expressed is None or n_features_expressed <= 0:
		return np.nan
	ratio = float(anticor_gene_count) / float(n_features_expressed)
	ratio = max(0.0, ratio)
	return float(np.exp(-float(score_k) * ratio))


def run_anticor_per_celltype(
	h5ad_path,
	output_csv,
	cluster_key=None,
	min_cells=10,
	max_cell_types=None,
	species="hsapiens",
	n_rand_feat=2000,
	fpr=0.001,
	fdr=1 / 15,
	num_pos_cor=10,
	bin_size=5000,
	scratch_dir=None,
	offline_mode=False,
	use_live_pathway_lookup=False,
	score_k=1.0,
):
	"""Run anticor feature selection per cell type and export one CSV."""
	adata = sc.read_h5ad(h5ad_path)
	cluster_col = _pick_cluster_key(adata, cluster_key)
	feature_ids = list(map(str, adata.var_names.tolist()))

	cluster_series = adata.obs[cluster_col].astype(str)
	cell_types = sorted(cluster_series.unique().tolist())
	if max_cell_types is not None:
		cell_types = cell_types[: int(max_cell_types)]

	rows = []
	input_path = str(Path(h5ad_path).resolve())

	for cell_type in cell_types:
		idx = np.where(cluster_series.values == cell_type)[0]
		n_cells = int(len(idx))

		row = {
			"input_h5ad": input_path,
			"cluster_key": cluster_col,
			"celltype": cell_type,
			"score": np.nan,
			"cell_type": cell_type,
			"n_cells": n_cells,
			"n_features_total": np.nan,
			"anticor_gene_count": np.nan,
			"anticor_score_0_1": np.nan,
			"status": "ok",
			"error": "",
		}

		if n_cells < int(min_cells):
			row["status"] = "skipped"
			row["error"] = f"n_cells < min_cells ({n_cells} < {min_cells})"
			rows.append(row)
			continue

		try:
			x_sub = adata.X[idx, :]
			expressed_mask = _expressed_gene_mask(x_sub)
			expressed_genes_set = set(np.asarray(feature_ids, dtype=object)[expressed_mask].tolist())
			n_features_expressed = _count_expressed_genes(x_sub)
			row["n_features_total"] = int(n_features_expressed)
			exprs_gene_by_cell = _to_gene_by_cell_matrix(x_sub)

			feature_table = anticor_pkg.get_anti_cor_genes(
				exprs=exprs_gene_by_cell,
				feature_ids=feature_ids,
				species=species,
				n_rand_feat=int(n_rand_feat),
				FPR=float(fpr),
				FDR=float(fdr),
				num_pos_cor=int(num_pos_cor),
				bin_size=int(bin_size),
				scratch_dir=scratch_dir,
				cell_axis=1,
				offline_mode=bool(offline_mode),
				use_live_pathway_lookup=bool(use_live_pathway_lookup),
			)

			if "selected" in feature_table.columns:
				selected_mask = feature_table["selected"].astype(bool).to_numpy()
				selected_genes = feature_table.loc[selected_mask, "gene"].astype(str)
				n_selected = int(np.sum(selected_genes.isin(expressed_genes_set)))
			else:
				# Conservative fallback if package output schema changes.
				selected_genes = feature_table["gene"].astype(str)
				n_selected = int(np.sum(selected_genes.isin(expressed_genes_set)))

			row["anticor_gene_count"] = n_selected
			row["anticor_score_0_1"] = _asymptotic_score(
				anticor_gene_count=n_selected,
				n_features_expressed=n_features_expressed,
				score_k=score_k,
			)
			row["score"] = row["anticor_score_0_1"]
		except Exception as exc:
			row["status"] = "failed"
			row["error"] = str(exc)

		rows.append(row)

	out_df = pd.DataFrame(rows)

	out_df = out_df.sort_values(
		["anticor_score_0_1", "anticor_gene_count", "n_cells"],
		ascending=[False, False, False],
		na_position="last",
	)

	output_path = Path(output_csv)
	output_path.parent.mkdir(parents=True, exist_ok=True)
	out_df.to_csv(output_path, index=False)
	return output_path


def main():
	parser = argparse.ArgumentParser(
		description="Run anti-correlation feature discovery per cell type and export CSV."
	)
	parser.add_argument("--input", required=True, help="Path to input .h5ad file")
	parser.add_argument(
		"--output",
		default=None,
		help="Output CSV path (default: <input_stem>_anticor_per_cluster.csv)",
	)
	parser.add_argument(
		"--cluster-key",
		default=None,
		help="obs column with cell type/cluster labels (auto-detect if omitted)",
	)
	parser.add_argument(
		"--min-cells",
		type=int,
		default=10,
		help="Minimum cells required for a cell type to run analysis (default: 10)",
	)
	parser.add_argument(
		"--max-cell-types",
		type=int,
		default=None,
		help="Optional cap on number of cell types (useful for quick testing)",
	)
	parser.add_argument(
		"--species",
		default="hsapiens",
		help="Species argument passed to get_anti_cor_genes (default: hsapiens)",
	)
	parser.add_argument(
		"--n-rand-feat",
		type=int,
		default=2000,
		help="n_rand_feat passed to get_anti_cor_genes (default: 2000)",
	)
	parser.add_argument("--fpr", type=float, default=0.001, help="FPR threshold")
	parser.add_argument("--fdr", type=float, default=1 / 15, help="FDR threshold")
	parser.add_argument(
		"--num-pos-cor",
		type=int,
		default=10,
		help="num_pos_cor passed to get_anti_cor_genes",
	)
	parser.add_argument(
		"--bin-size",
		type=int,
		default=5000,
		help="bin_size passed to get_anti_cor_genes",
	)
	parser.add_argument(
		"--scratch-dir",
		default=None,
		help="Optional scratch directory for package temp files",
	)
	parser.add_argument(
		"--offline-mode",
		action="store_true",
		help="Run with offline_mode=True in get_anti_cor_genes",
	)
	parser.add_argument(
		"--use-live-pathway-lookup",
		action="store_true",
		help="Enable live pathway lookup in get_anti_cor_genes",
	)
	parser.add_argument(
		"--score-k",
		type=float,
		default=1.0,
		help="Asymptotic score steepness k in exp(-k * ratio), where ratio=anticor/n_features_total",
	)

	args = parser.parse_args()
	input_path = Path(args.input)
	output_csv = args.output or str(
		input_path.with_name(f"{input_path.stem}_anticor_per_cluster.csv")
	)

	output_path = run_anticor_per_celltype(
		h5ad_path=str(input_path),
		output_csv=output_csv,
		cluster_key=args.cluster_key,
		min_cells=args.min_cells,
		max_cell_types=args.max_cell_types,
		species=args.species,
		n_rand_feat=args.n_rand_feat,
		fpr=args.fpr,
		fdr=args.fdr,
		num_pos_cor=args.num_pos_cor,
		bin_size=args.bin_size,
		scratch_dir=args.scratch_dir,
		offline_mode=args.offline_mode,
		use_live_pathway_lookup=args.use_live_pathway_lookup,
		score_k=args.score_k,
	)
	print(f"Saved anti-correlation per-cell-type scores to: {output_path}")


if __name__ == "__main__":
	main()
