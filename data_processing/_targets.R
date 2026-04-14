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
  # Track config files so changes invalidate downstream targets
  tar_target(
    processing_parameters_file,
    "config/processing_parameters.yaml",
    format = "file"
  ),
  tar_target(
    core_datasets_file,
    "config/core_datasets.yaml",
    format = "file"
  ),

  # Load configuration once
  tar_target(
    config,
    {
      processing_parameters_file
      load_dp_params()
    }
  ),
  
  # Load dataset registry once
  tar_target(
    datasets,
    {
      core_datasets_file
      load_dp_datasets()
    }
  ),
  
  # Dynamically branch once per dataset in the registry.
  # This produces branch targets like process_ds_1, process_ds_2, ...
  tar_target(
    process_ds,
    {
      # In this targets version, each branch receives a single-row object
      # named `datasets` (not individual column variables).
      row <- datasets
      id <- row$id
      reference <- row$reference
      raw_file <- row$raw_file[[1]]
      raw_files <- row$raw_files[[1]]
      ident_col <- row$ident_col
      sample_col <- row$sample_col
      batch_col <- row$batch_col
      condition_col <- row$condition_col
      exclude_cell_types <- row$exclude_cell_types[[1]]

      # Reconstitute a dataset config list for helper functions.
      ds <- list(
        id = id,
        reference = reference,
        raw_file = raw_file,
        raw_files = raw_files,
        ident_col = ident_col,
        sample_col = sample_col,
        batch_col = batch_col,
        condition_col = if (!is.null(condition_col) && nzchar(condition_col)) condition_col else NULL,
        exclude_cell_types = exclude_cell_types
      )

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
      obj <- load_dataset(ds, input_dir = config$in_dir)

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
    pattern = map(datasets),
    iteration = "list",
    # CRITICAL: Continue on error so other datasets aren't affected
    error = "continue"
  ),
  
  # Summary target: collect results from all datasets
  tar_target(
    summary,
    {
      # `process_ds` is a list of per-branch results.
      results <- do.call(
        rbind,
        lapply(process_ds, function(x) as.data.frame(x, stringsAsFactors = FALSE))
      )
      rownames(results) <- NULL
      results
    }
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
