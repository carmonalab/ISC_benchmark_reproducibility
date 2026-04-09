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

# Source utilities from parent directory and local modules
source("../R/cli_utils.R")
source("../R/shared_helpers.R")
source("R/00_utils.R")
source("R/data_processing_helpers.R")
source("R/process_datasets.R")

# ============================================================================
# TARGETS PIPELINE CONFIGURATION
# ============================================================================

tar_option_set(
  packages = c("yaml", "dplyr", "stringr", "Seurat", "SeuratObject", "Matrix", "BiocParallel"),
  # Store outputs in data/processed/
  storage = "worker",
  retrieval = "worker"
)

# ============================================================================
# DYNAMIC TARGETS: ONE PER DATASET
# ============================================================================

list(
  # Load configuration once
  tar_target(
    config,
    load_dp_params(),
    # Invalidate if config file changes
    cue = tar_cue(file = "config/processing_parameters.yaml")
  ),
  
  # Load dataset registry once
  tar_target(
    datasets,
    load_dpd_dp_datasets(),
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
