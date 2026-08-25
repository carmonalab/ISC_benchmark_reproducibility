#!/usr/bin/env python3
"""Run SCCAF on an h5ad file and export per-cluster metrics as CSV.

Workflow:
1. Load AnnData from .h5ad.
2. Pick a clustering column from adata.obs (manual or auto-detected).
3. Run SCCAF_assessment once.
4. Compute per-cluster precision, recall, and F1 from SCCAF hold-out predictions.
5. Save one CSV with per-cluster values.
"""

import argparse
import csv
from collections import Counter
from pathlib import Path

import scanpy as sc
from SCCAF import SCCAF_assessment
from sklearn.metrics import precision_recall_fscore_support


def _pick_cluster_key(adata, requested_key=None):
	"""Return cluster-label column name from adata.obs."""
	if requested_key:
		if requested_key not in adata.obs.columns:
			raise ValueError(
				f"Cluster column '{requested_key}' not found in adata.obs. "
				f"Available columns: {list(adata.obs.columns)}"
			)
		return requested_key

	preferred = [
		"louvain",
		"leiden",
		"seurat_clusters",
		"OriginalAnnotationLevel2",
		"OriginalAnnotationLevel1",
		"cell_type",
		"celltype",
		"annotation",
	]
	for key in preferred:
		if key in adata.obs.columns:
			return key

	raise ValueError(
		"No cluster column provided and none of the default columns were found. "
		f"Available columns: {list(adata.obs.columns)}"
	)


def run_sccaf_per_cluster(h5ad_path, output_csv, cluster_key=None, n=100):
	"""Run SCCAF and export per-cluster metrics.

	Parameters
	----------
	h5ad_path : str
		Path to input h5ad file.
	output_csv : str
		Path to per-cluster output CSV.
	cluster_key : str or None
		Column in adata.obs with cluster labels; if None, auto-detect.
	n : int
		SCCAF assessment iterations (passed to SCCAF_assessment).
	"""
	adata = sc.read_h5ad(h5ad_path)
	key = _pick_cluster_key(adata, cluster_key)

	labels = adata.obs[key].astype(str)
	_, y_pred, y_test, _, _, _ = SCCAF_assessment(adata.X, labels, n=n)

	output_path = Path(output_csv)
	output_path.parent.mkdir(parents=True, exist_ok=True)

	true_labels = [str(x) for x in y_test]
	pred_labels = [str(x) for x in y_pred]

	if len(true_labels) != len(pred_labels):
		raise ValueError(
			"SCCAF output mismatch: y_test and y_pred have different lengths "
			f"({len(true_labels)} vs {len(pred_labels)})."
		)

	# input counts use all cells before SCCAF hold-out split
	input_counts = Counter([str(x) for x in labels.values])
	true_counts = Counter(true_labels)
	pred_counts = Counter(pred_labels)
	correct_counts = Counter(
		true for true, pred in zip(true_labels, pred_labels) if true == pred
	)

	all_clusters = sorted(set(true_labels) | set(pred_labels))
	precision, recall, f1, support = precision_recall_fscore_support(
		true_labels,
		pred_labels,
		labels=all_clusters,
		zero_division=0,
	)

	with output_path.open("w", newline="") as handle:
		writer = csv.DictWriter(
			handle,
			fieldnames=[
				"input_h5ad",
				"cluster_key",
				"cluster",
				"n_cells_input",
				"n_cells_test",
				"n_predicted_as_cluster",
				"n_correct",
				"recall_within_cluster",
				"precision_for_cluster",
				"f1_score",
			],
		)
		writer.writeheader()

		for i, cluster in enumerate(all_clusters):
			test_support = int(support[i])
			predicted = int(pred_counts.get(cluster, 0))
			correct = int(correct_counts.get(cluster, 0))

			writer.writerow(
				{
					"input_h5ad": str(Path(h5ad_path).resolve()),
					"cluster_key": key,
					"cluster": cluster,
					"n_cells_input": int(input_counts.get(cluster, 0)),
					"n_cells_test": test_support,
					"n_predicted_as_cluster": predicted,
					"n_correct": correct,
					"recall_within_cluster": float(recall[i]),
					"precision_for_cluster": float(precision[i]),
					"f1_score": float(f1[i]),
				}
			)

	return output_path


def main():
	parser = argparse.ArgumentParser(
		description="Run SCCAF on an h5ad and export per-cluster metrics CSV."
	)
	parser.add_argument(
		"--input",
		required=True,
		help="Path to input .h5ad file",
	)
	parser.add_argument(
		"--output",
		default=None,
		help="Path to per-cluster output CSV (default: <input_stem>_sccaf_per_cluster.csv)",
	)
	parser.add_argument(
		"--cluster-key",
		default=None,
		help="Column in adata.obs with cluster labels (default: auto-detect)",
	)
	parser.add_argument(
		"--n",
		type=int,
		default=100,
		help="Number of SCCAF assessment iterations (default: 100)",
	)

	args = parser.parse_args()

	input_path = Path(args.input)
	output_csv = args.output or str(input_path.with_name(f"{input_path.stem}_sccaf_per_cluster.csv"))

	out_path = run_sccaf_per_cluster(
		h5ad_path=str(input_path),
		output_csv=output_csv,
		cluster_key=args.cluster_key,
		n=args.n,
	)
	print(f"Saved SCCAF per-cluster metrics to: {out_path}")


if __name__ == "__main__":
	main()
