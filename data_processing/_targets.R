# data_processing/_targets.R
#
# Targets pipeline for dataset processing
#
# Features:
#   - Incremental processing: only reruns if outputs are missing
#   - Error handling: if one dataset fails, others continue
#   - Dependency tracking: automatically detects config/source changes
#   - Progress visibility: shows which datasets succeeded/failed
#
# Usage:
#   cd data_processing
#   targets::tar_make()                    # Run all pending targets
#   targets::tar_make(asks = FALSE)        # Force rerun all
#   targets::tar_cue(names = "config")     # Rerun if config changed
#   targets::tar_status()                  # Check status
#
# NOTE: Run from data_processing directory (NOT from project root)

library(targets)
library(yaml)

# Source utilities from parent directory
source("../R/project_paths.R")
source("../R/cli_utils.R")
source("../R/shared_helpers.R")
source("R/data_processing_helpers.R")

# ============================================================================
# TARGETS PIPELINE CONFIGURATION
# ============================================================================

tar_option_set(
  packages = c("yaml", "dplyr", "stringr", "Seurat", "SeuratObject", "Matrix", "BiocParallel"),
  # Store outputs in data/processed/isc/
  storage = "worker",
  retrieval = "worker"
)

# ============================================================================
# HELPER FUNCTIONS FOR TARGETS
# ============================================================================

#' Load configuration
#' @return List of configuration parameters
get_config <- function() {
  config <- read_yaml("config/processing_parameters.yaml")
  config$out_dir <- "../data/processed/isc"  # Relative to data_processing dir
  dir.create(config$out_dir, showWarnings = FALSE, recursive = TRUE)
  config
}

#' Load dataset registry
#' @return List of dataset configurations
get_datasets <- function() {
  datasets_config <- read_yaml("config/core_datasets.yaml")
  datasets_config$datasets
}

#' Check if dataset output files exist
#' @param prefix Dataset prefix (e.g., "JoaI")
#' @param out_dir Output directory
#' @return TRUE if at least one output file exists, FALSE otherwise
output_exists <- function(prefix, out_dir) {
  files <- list.files(out_dir, pattern = paste0("^", prefix, "_.*\\.rds$"))
  length(files) > 0
}

# ============================================================================
# DYNAMIC TARGETS: ONE PER DATASET
# ============================================================================

list(
  # Load configuration once
  tar_target(
    config,
    get_config(),
    # Invalidate if config file changes
    cue = tar_cue(file = "config/processing_parameters.yaml")
  ),
  
  # Load dataset registry once
  tar_target(
    datasets,
    get_datasets(),
    # Invalidate if registry changes
    cue = tar_cue(file = "config/core_datasets.yaml")
  ),
  
  # Dynamically create targets for each dataset
  tar_map(
    # Iterate over datasets
    values = datasets,
    
    # Per-dataset target
    tar_target(
      name = process_ds,
      command = {
        # Extract dataset ID (remove year suffix)
        prefix <- sub("_[0-9]+$", "", id)
        
        # Check if output already exists
        if (output_exists(prefix, config$out_dir)) {
          message(sprintf("[SKIP] %s - outputs already exist", reference))
          return(list(
            prefix = prefix,
            status = "skipped",
            reason = "outputs_exist",
            n_files = length(list.files(config$out_dir, pattern = paste0("^", prefix, "_.*\\.rds$")))
          ))
        }
        
        message(sprintf("[PROCESSING] %s", reference))
        
        # Load dataset
        obj <- load_dataset(
          list(
            id = id,
            raw_file = raw_file,
            raw_files = raw_files,
            reference = reference
          )
        )
        
        if (is.null(obj)) {
          message(sprintf("  ERROR: Failed to load dataset %s", id))
          return(list(
            prefix = prefix,
            status = "failed",
            reason = "load_failed",
            n_files = 0
          ))
        }
        
        # Apply dataset-specific preprocessing
        obj <- tryCatch({
          # Create dataset config object from tar_map values
          ds <- list(
            id = id,
            raw_file = raw_file,
            raw_files = raw_files,
            reference = reference
          )
          preprocess_dataset(obj, ds)
        }, error = function(e) {
          message(sprintf("  ERROR during preprocessing: %s", e$message))
          return(NULL)
        })
        
        if (is.null(obj)) {
          return(list(
            prefix = prefix,
            status = "failed",
            reason = "preprocessing_failed",
            n_files = 0
          ))
        }
        
        # Process dataset
        result <- tryCatch({
          process_dataset(
            obj = obj,
            prefix = prefix,
            ident_col = ident_col,
            sample_col = sample_col,
            batch_col = batch_col,
            condition_col = condition_col,
            exclude_celltypes = exclude_cell_types,
            config = config
          )
          
          # Count output files
          n_files <- length(list.files(config$out_dir, pattern = paste0("^", prefix, "_.*\\.rds$")))
          
          list(
            prefix = prefix,
            status = "success",
            reason = NA,
            n_files = n_files
          )
        }, error = function(e) {
          message(sprintf("  ERROR during processing: %s", e$message))
          list(
            prefix = prefix,
            status = "failed",
            reason = "processing_failed",
            n_files = 0
          )
        })
        
        rm(obj)
        gc()
        
        result
      },
      # CRITICAL: Continue on error so other datasets aren't affected
      error = "continue"
    )
  ),
  
  # Summary target: collect results from all datasets
  tar_target(
    summary,
    {
      # Bind all dataset results into one data frame
      results <- do.call(rbind, lapply(
        grep("^process_ds_", names(.targets), value = TRUE),
        get
      ))
      
      rownames(results) <- NULL
      results
    },
    # Depends on all dataset processing targets
    pattern = cross(grep("^process_ds_", names(.targets), value = TRUE))
  ),
  
  # Report target
  tar_target(
    report,
    {
      cat("\n")
      cat("="*60, "\n")
      cat("DATA PROCESSING SUMMARY\n")
      cat("="*60, "\n")
      cat("\n")
      
      # Summary table
      cat("Status Summary:\n")
      print(table(summary$status))
      cat("\n")
      
      # Detailed results
      cat("Detailed Results:\n")
      for (i in seq_len(nrow(summary))) {
        row <- summary[i, ]
        status_icon <- if (row$status == "success") "✓" else if (row$status == "skipped") "-" else "✗"
        cat(sprintf("%s %-20s: %s (outputs: %d)\n",
                   status_icon,
                   row$prefix,
                   row$status,
                   row$n_files))
        if (!is.na(row$reason)) {
          cat(sprintf("  └─ %s\n", row$reason))
        }
      }
      
      cat("\n")
      cat("="*60, "\n")
      
      # Total counts
      success_count <- sum(summary$status == "success")
      failed_count <- sum(summary$status == "failed")
      skipped_count <- sum(summary$status == "skipped")
      total_files <- sum(summary$n_files)
      
      cat(sprintf("Successful: %d | Failed: %d | Skipped: %d\n", 
                 success_count, failed_count, skipped_count))
      cat(sprintf("Total output files: %d\n", total_files))
      cat("="*60, "\n\n")
      
      invisible(summary)
    }
  )
)
