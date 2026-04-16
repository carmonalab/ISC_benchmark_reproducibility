#!/usr/bin/env Rscript

#' Prepare Label Transfer Query/Reference Splits
#'
#' Post-processing script to create query (70% samples) and reference (30% samples) 
#' splits from ISC-processed Seurat objects, following:
#' 
#'   Consistency_metrics_benchmark/OPS/manual_classifiers/process_datasets.R
#'
#' Key features:
#'   - Loads Seurat objects from data/processed/
#'   - Extracts counts and metadata using scTypeEval::load_singleCell_object
#'   - Stratified sample-level split (seed from label_transfer_parameters.yaml)
#'   - Supports multiple split replicates (n_replicates) with different seeds
#'   - Skips datasets with <min_samples samples (from config)
#'   - Outputs query.rds and reference.rds to data/processed/label_transfer/rep<k>/<dataset_id>/
#'
#' Usage:
#'   Rscript label_transfer_task/R/00_prepare_splits.R
#'
#' Requirements:
#'   - Run from project root: cd /path/to/ISC_benchmark_reproducibility
#'   - ISC processed datasets already in data/processed/
#'   - R packages: dplyr, scTypeEval, Seurat, BiocParallel
#'

suppressPackageStartupMessages({
  library(BiocParallel)
  library(yaml)
})

# Source shared utilities (robust to being run from any working directory)
get_script_path <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  if (length(file_arg) > 0) {
    return(sub("^--file=", "", file_arg[[1]]))
  }
  NULL
}

script_path <- get_script_path()
if (is.null(script_path)) {
  stop("Unable to determine script path (expected Rscript --file=...)")
}

script_dir <- dirname(normalizePath(script_path))
proj_root_guess <- normalizePath(file.path(script_dir, "../.."))

source(file.path(proj_root_guess, "R", "shared_helpers.R"))
source(file.path(proj_root_guess, "R", "cli_utils.R"))

# Ensure the pipeline context is label_transfer_task for config loading
setwd(file.path(proj_root_guess, "label_transfer_task"))

source("R/00_utils.R")

# ============================================================================
# Configuration (from label_transfer_parameters.yaml)
# ============================================================================

params <- load_lt_params()
n_replicates <- get_lt_n_replicates(params)
n_cores <- get_lt_n_cores(params)

main <- function() {
  # Verify we're in correct directory
  if (!file.exists(proj_path("data_processing/config/datasets.yaml")) ||
      !dir.exists(proj_path("data/processed"))) {
    stop("Project data not found. Expected to find data_processing/config/datasets.yaml and data/processed under the project root.")
  }
  
  message_time(strrep("=", 70))
  message_time("Preparing Label Transfer Query/Reference Splits")
  message_time(strrep("=", 70))
  
  dataset_ids <- list_label_transfer_datasets_from_isc(params)
  if (length(dataset_ids) == 0) {
    stop("No eligible ISC datasets found for label transfer")
  }

  message_time(paste("Preparing splits for", length(dataset_ids), "datasets and", n_replicates, "replicates..."))

  workers <- min(n_cores, length(dataset_ids))

  all_results <- list()
  for (rep in seq_len(n_replicates)) {
    message_time(paste("--- Replicate", rep, "---"))
    res <- bplapply(
      dataset_ids,
      BPPARAM = BiocParallel::MulticoreParam(workers = workers, progressbar = TRUE),
      function(dataset_id) {
        tryCatch(
          prepare_label_transfer_split(dataset_id = dataset_id, replicate = rep, params = params),
          error = function(e) {
            warning("Error preparing ", dataset_id, " (rep ", rep, "): ", conditionMessage(e))
            NULL
          }
        )
      }
    )
    all_results[[paste0("rep", rep)]] <- Filter(Negate(is.null), res)
  }
  
  message_time(strrep("=", 70))
  message_time("✓ Split preparation complete")
  message_time(strrep("=", 70))
  message_time("")
  message_time("Label transfer splits saved to: data/processed/label_transfer/rep<k>/")
  message_time("")
  message_time("Next steps:")
  message_time("  1. Review splits: ls -lh data/processed/label_transfer/")
  message_time("  2. Run label transfer benchmark:")
  message_time("     cd label_transfer_task && Rscript -e 'targets::tar_make()'")
  message_time("")
  
  invisible(all_results)
}

if (!interactive()) {
  main()
}
