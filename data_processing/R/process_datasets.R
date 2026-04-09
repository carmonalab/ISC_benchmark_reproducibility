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
  
  results <- list()
  
  for (i in seq_along(datasets)) {
    ds <- datasets[[i]]
    prefix <- sub("_[0-9]+$", "", ds$id)  # Remove year suffix for cleaner names
    
    message(sprintf("Dataset %d/%d: %s", i, length(datasets), ds$reference))
    
    # Wrap entire dataset processing in try-catch for error resilience
    result <- tryCatch({
      # Load dataset
      obj <- load_dataset(ds)
      if (is.null(obj)) {
        message(sprintf("  SKIPPED (file not found)\n"))
        return(list(
          prefix = prefix,
          status = "skipped",
          n_files = 0,
          error = "file_not_found"
        ))
      }
      
      # Apply preprocessing (foreseen error: dataset-specific issues)
      obj <- tryCatch({
        preprocess_dataset(obj, ds)
      }, error = function(e) {
        stop(sprintf("Preprocessing failed: %s", e$message))
      })
      
      # Process dataset (foreseen error: filtering, downsampling, etc.)
      tryCatch({
        process_dataset(
          obj = obj,
          prefix = prefix,
          ident_col = ds$ident_col,
          sample_col = ds$sample_col,
          batch_col = ds$batch_col,
          condition_col = ds$condition_col,
          exclude_celltypes = ds$exclude_cell_types,
          config = config
        )
      }, error = function(e) {
        stop(sprintf("Processing failed: %s", e$message))
      })
      
      # Count output files
      n_files <- length(list.files(config$out_dir, pattern = paste0("^", prefix, "_.*\\.rds$")))
      
      rm(obj)
      gc()
      
      list(
        prefix = prefix,
        status = "success",
        n_files = n_files,
        error = NA
      )
    }, error = function(e) {
      message(sprintf("  ERROR: %s\n", e$message))
      list(
        prefix = prefix,
        status = "failed",
        n_files = 0,
        error = e$message
      )
    })
    
    results[[i]] <- result
    message("")
  }
  
  # Summary
  message("")
  message(sprintf("[%s] Processing Complete", format(Sys.time(), "%H:%M:%S")))
  message(sprintf("Output directory: %s", config$out_dir))
  message("")
  
  # Summarize results
  results_df <- do.call(rbind, lapply(results, as.data.frame))
  
  cat("==========================================================\n")
  cat("PROCESSING SUMMARY\n")
  cat("==========================================================\n")
  cat("\n")
  
  # Count by status
  status_table <- table(results_df$status)
  cat("Status counts:\n")
  print(status_table)
  cat("\n")
  
  # Detailed results
  cat("Detailed results:\n")
  for (i in seq_len(nrow(results_df))) {
    row <- results_df[i, ]
    status_icon <- if (row$status == "success") "✓" else if (row$status == "skipped") "-" else "✗"
    cat(sprintf("%s %-20s: %s (outputs: %d)\n",
               status_icon,
               row$prefix,
               row$status,
               row$n_files))
    if (!is.na(row$error) && row$status == "failed") {
      cat(sprintf("  └─ Error: %s\n", row$error))
    }
  }
  
  cat("\n")
  cat(sprintf("Total output files: %d\n", sum(results_df$n_files)))
  cat("==========================================================\n\n")
  
  invisible(results_df)
}

# ============================================================================
# Entry Point
# ============================================================================

if (!interactive()) {
  main()
}
