#!/usr/bin/env Rscript

#' ISC Benchmark Data Processing Pipeline
#'
#' Reproduces the data processing from Consistency_metrics_benchmark/datasets/proc_data.Rmd
#' Processes 6 core datasets, downsamples cells, splits by batch + condition,
#' and applies optimal_dataset filtering.
#'
#' Configuration: data_processing/config/processing_parameters.yaml
#' Dataset registry: data_processing/config/core_datasets.yaml
#'
#' Output structure: data/processed/isc/PREFIX_batch_condition.rds
#'
#' Usage:
#'   Rscript data_processing/process_datasets.R [--dry-run]

suppressPackageStartupMessages({
  library(yaml)
  library(dplyr)
  library(stringr)
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(BiocParallel)
})

# Load utilities
source("R/project_paths.R")
source("R/cli_utils.R")
source("R/shared_helpers.R")
source("data_processing/R/data_processing_helpers.R")

# ============================================================================
# LOAD CONFIGURATION FROM YAML
# ============================================================================

config <- read_yaml("data_processing/config/processing_parameters.yaml")

# Set output directory
dir.create(config$out_dir, showWarnings = FALSE, recursive = TRUE)

message(sprintf("\n[%s] ISC Benchmark Data Processing", format(Sys.time(), "%H:%M:%S")))
message(sprintf("  Seed: %d", config$seed))
message(sprintf("  Cores: %d", config$n_cores))
message(sprintf("  Output: %s", config$out_dir))
message("")

# ============================================================================
# LOAD DATASET REGISTRY FROM YAML
# ============================================================================

datasets_config <- read_yaml("data_processing/config/core_datasets.yaml")
datasets <- datasets_config$datasets

message(sprintf("Loaded %d datasets from core_datasets.yaml", length(datasets)))
message("")

# ============================================================================
# MAIN: ITERATE THROUGH DATASETS
# ============================================================================

main <- function() {
  message(sprintf("Beginning processing of %d datasets...\n", length(datasets)))
  
  for (i in seq_along(datasets)) {
    ds <- datasets[[i]]
    
    message(sprintf("Dataset %d/%d: %s", i, length(datasets), ds$reference))
    
    # Load dataset
    obj <- load_dataset(ds)
    if (is.null(obj)) {
      message(sprintf("  SKIPPED (file not found)\n"))
      next
    }
    
    # Apply preprocessing
    obj <- preprocess_dataset(obj, ds)
    
    # Process dataset
    process_dataset(
      obj = obj,
      prefix = sub("_[0-9]+$", "", ds$id),  # Remove year suffix for cleaner names
      ident_col = ds$ident_col,
      sample_col = ds$sample_col,
      batch_col = ds$batch_col,
      condition_col = ds$condition_col,
      exclude_celltypes = ds$exclude_cell_types,
      config = config
    )
    
    rm(obj)
    gc()
    message("")
  }
  
  # Summary
  message("")
  message(sprintf("[%s] Processing Complete", format(Sys.time(), "%H:%M:%S")))
  message(sprintf("Output directory: %s", config$out_dir))
  
  # Count output files
  out_files <- list.files(config$out_dir, pattern = "\\.rds$")
  message(sprintf("Total output files: %d", length(out_files)))
  
  if (length(out_files) > 0) {
    message("\nOutput files:")
    invisible(lapply(head(out_files, 10), function(f) message(sprintf("  - %s", f))))
    if (length(out_files) > 10) {
      message(sprintf("  ... and %d more", length(out_files) - 10))
    }
  }
  
  invisible(NULL)
}

# ============================================================================
# Entry Point
# ============================================================================

if (!interactive()) {
  main()
}
