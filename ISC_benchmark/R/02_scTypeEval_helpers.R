# ISC_benchmark/R/02_scTypeEval_helpers.R --- scTypeEval metric integration
#
# Wraps scTypeEval functions for ISC benchmark metric computation.
# Uses utilities from scTypeEval/inst/benchmarking/{Metrics_benchmarking.R, assays_utils.R}
#
# Integration points with scTypeEval's wr_* wrappers:
# ├─ wr_missclasify:     Computes metrics on shuffled label perturbations (Task 1)
# ├─ wr_nsamples:        Computes metrics on downsampled datasets (Task 2)
# ├─ wr_nct:            Computes metrics on cell-type combinations (Task 3)
# ├─ wr_ncell:          Computes metrics on cell-downsampled datasets (Task 4)
# ├─ wr_merge_ct:       Computes metrics on merged cell type annotations (Task 5)
# └─ wr_split_cell_type: Computes metrics on split cell type annotations (Task 6)
#
# wr_* functions use common patterns:
#  1. wrapper_dissimilarity() — Runs full dissimilarity pipeline
#  2. get_consistency() — Extracts consistency metrics from scTypeEval object
#  3. Compares original_ident vs perturbed annotations

# ============================================================================
# LOAD scTypeEval BENCHMARKING UTILITIES
# ============================================================================

load_scTypeEval_benchmarking <- function() {
  # Try to load scTypeEval's benchmarking utilities
  
  if (!requireNamespace("scTypeEval", quietly = TRUE)) {
    stop("scTypeEval not found. Install via: ",
         "remotes::install_github('carmonalab/scTypeEval@v0.99.30')")
  }
  
  # Path to scTypeEval's built-in benchmarking scripts
  metrics_file <- tryCatch(
    system.file("benchmarking", "Metrics_benchmarking.R", package = "scTypeEval"),
    error = function(e) ""
  )
  
  assays_file <- tryCatch(
    system.file("benchmarking", "assays_utils.R", package = "scTypeEval"),
    error = function(e) ""
  )
  
  if (metrics_file == "" || assays_file == "") {
    warning("Could not locate scTypeEval benchmarking utilities. ",
            "Check installation: scTypeEval/inst/benchmarking/")
    return(FALSE)
  }
  
  # Source them into current environment
  try(source(assays_file, local = TRUE), silent = TRUE)
  try(source(metrics_file, local = TRUE), silent = TRUE)
  
  TRUE
}

load_scTypeEval_benchmarking()
