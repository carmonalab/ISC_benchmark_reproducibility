#!/usr/bin/env Rscript
# ==============================================================================
# Run All Cell Type Classifiers
# ==============================================================================
# Main script that applies all classifiers to query datasets
# Features:
#   - Loads classifier functions from classifiers.R
#   - Runs classifiers with error handling (continue if one fails)
#   - Incremental execution (skips already-completed datasets)
#   - Saves predictions.rds with classifier predictions per cell
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(BiocParallel)
  library(yaml)
})

# Load utilities
# First establish project root using shared_helpers function
source("../R/shared_helpers.R")
source(proj_path("R/cli_utils.R"))  # For message functions
source("R/00_utils.R")              # For load_lt_params, get_idents_by_prefix, etc.

# ============================================================================
# CONFIGURATION
# ============================================================================

params <- load_lt_params()
SEED <- get_lt_seed(params)
N_CORES <- get_lt_n_cores(params)

# ============================================================================
# HELPER: Check if dataset already completed
# ============================================================================

#' Check if classifier predictions already exist for dataset
#'
#' @param dataset_dir Path to dataset directory (contains query.rds, reference.rds)
#' @param strict If TRUE, require all classifiers completed; if FALSE, any is fine
#'
#' @return TRUE if predictions.rds file exists, FALSE otherwise
is_dataset_completed <- function(dataset_dir, strict = FALSE) {
  predictions_file <- file.path(dataset_dir, "predictions.rds")
  file.exists(predictions_file)
}

# ============================================================================
# CLASSIFIER EXECUTION
# ============================================================================

#' Run all classifiers on a single dataset
#'
#' @param dataset_dir Directory containing query.rds and reference.rds
#' @param dataset_name Name of dataset (for logging)
#' @param params Pipeline parameters
#'
#' @return Updated metadata with classifier predictions or NULL if failed
classify_dataset <- function(dataset_dir, dataset_name, params, ncores = 1) {
  
  message_step("CLASSIFY", sprintf("Processing dataset: %s", dataset_name))
  
  # Check if already completed
  if (is_dataset_completed(dataset_dir)) {
    message_step("SKIP", sprintf("Predictions already exist: %s", dataset_name))
    return(readRDS(file.path(dataset_dir, "predictions.rds")))
  }
  
  # Load data
  query_file <- file.path(dataset_dir, "query.rds")
  reference_file <- file.path(dataset_dir, "reference.rds")
  
  if (!file.exists(query_file) || !file.exists(reference_file)) {
    message_step("ERROR", sprintf("Missing query or reference files: %s", dataset_name))
    return(NULL)
  }
  
  query <- readRDS(query_file)
  reference <- readRDS(reference_file)
  
  ref_counts <- reference$counts
  ref_labels <- reference$metadata$cell_type
  query_counts <- query$counts
  query_md <- query$metadata
  
  message_step("INFO", sprintf(
    "Reference: %d cells (%d cell types) | Query: %d cells",
    ncol(ref_counts), n_distinct(ref_labels), ncol(query_counts)
  ))
  
  # Initialize predictions dataframe
  predictions_list <- list()
  
  # Get list of classifiers to run
  classifiers <- get_lt_classifiers(params)
  
  # Load classifier functions
  ncores <- ncores  # make available to classifiers.R
  source(proj_path("label_transfer_task/classifiers/classifiers.R"), local = environment())
  
  # Run each classifier
  for (clf_name in classifiers) {
    clf_func_name <- paste0("classify_", clf_name)
    
    # Check if classifier function exists
    if (!exists(clf_func_name, envir = environment())) {
      message_step("SKIP", sprintf("[%s] NOT FOUND", clf_name))
      next
    }
    
    message(sprintf("  [%s] ", clf_name), appendLF = FALSE)
    
    # Try to run classifier
    pred <- tryCatch({
      clf_func <- get(clf_func_name, envir = environment())
      
      # Special handling for some classifiers
      result <- clf_func(ref_counts, ref_labels, query_counts)
      
      if (is.na(result)) {
        rep(NA, ncol(query_counts))
      } else {
        result
      }
    }, error = function(e) {
      message(sprintf("ERROR: %s", e$message))
      NA
    })
    
    # Validate prediction
    if (!is.na(pred[1])) {
      if (length(pred) == ncol(query_counts)) {
        # Calculate accuracy if ground truth available
        accuracy <- mean(pred == query_md$cell_type, na.rm = TRUE)
        n_correct <- sum(pred == query_md$cell_type, na.rm = TRUE)
        n_na <- sum(is.na(pred))
        message(sprintf(
          "✓ %.1f%% accuracy (%d/%d, %d NA)",
          accuracy * 100,
          n_correct,
          sum(!is.na(pred)),
          n_na
        ))
      } else {
        message(sprintf("ERROR: length mismatch (%d != %d)", length(pred), ncol(query_counts)))
        pred <- rep(NA, ncol(query_counts))
      }
    } else {
      message("FAILED")
      pred <- rep(NA, ncol(query_counts))
    }
    
    predictions_list[[clf_name]] <- pred
  }
  
  # ========================================================================
  # ENSEMBLE VOTING (exclude Random classifier)
  # ========================================================================
  
  message(sprintf("  [Ensemble] "), appendLF = FALSE)
  
  predictions_for_ensemble <- predictions_list[names(predictions_list) != "Random"]
  
  if (length(predictions_for_ensemble) > 0) {
    ensemble_pred <- tryCatch({
      # Ensemble: majority voting across classifiers
      pred_matrix <- do.call(rbind, predictions_for_ensemble)
      
      # Majority vote: most common prediction per cell
      ensemble <- apply(pred_matrix, 2, function(col) {
        if (all(is.na(col))) {
          return(NA)
        }
        # Get most frequent non-NA prediction
        mode_pred <- names(sort(table(col), decreasing = TRUE))[1]
        if (is.na(mode_pred)) NA else mode_pred
      })
      
      ensemble
    }, error = function(e) {
      message(sprintf("ERROR: %s", e$message))
      NA
    })
    
    if (!is.na(ensemble_pred[1])) {
      accuracy <- mean(ensemble_pred == query_md$cell_type, na.rm = TRUE)
      n_correct <- sum(ensemble_pred == query_md$cell_type, na.rm = TRUE)
      n_na <- sum(is.na(ensemble_pred))
      message(sprintf(
        "✓ %.1f%% accuracy (%d/%d, %d NA)",
        accuracy * 100,
        n_correct,
        sum(!is.na(ensemble_pred)),
        n_na
      ))
    } else {
      message("FAILED")
      ensemble_pred <- rep(NA, ncol(query_counts))
    }
    
    predictions_list$Ensemble <- ensemble_pred
  }
  
  # ========================================================================
  # SAVE RESULTS
  # ========================================================================
  
  # Add predictions to metadata
  updated_md <- query_md %>% as.data.frame()
  rownames(updated_md) <- rownames(query_md)
  
  for (clf_name in names(predictions_list)) {
    col_name <- paste0("pred_", clf_name)
    updated_md[[col_name]] <- predictions_list[[clf_name]]
  }
  
  # Save updated metadata with predictions
  output_file <- file.path(dataset_dir, "predictions.rds")
  saveRDS(updated_md, output_file)
  
  message_step("SAVE", sprintf("Predictions saved: %s", output_file))
  
  return(updated_md)
}

# ============================================================================
# BATCH PROCESSING
# ============================================================================

#' Run all classifiers on all datasets
#'
#' @param base_dir Directory containing all dataset subdirectories
#' @param params Pipeline parameters
#'
#' @return List of results per dataset
run_all_classifiers <- function(base_dir = NULL, params = NULL) {
  
  if (is.null(base_dir)) {
    base_dir <- proj_path("data/processed/label_transfer")
  }
  
  if (is.null(params)) {
    params <- load_lt_params()
  }
  
  base_dir <- normalizePath(base_dir, mustWork = FALSE)
  
  if (!dir.exists(base_dir)) {
    stop("Label-transfer data directory not found: ", base_dir)
  }
  
  # Support both layouts:
  #  - legacy: <base_dir>/<dataset_id>/{query,reference}.rds
  #  - replicate: <base_dir>/rep<k>/<dataset_id>/{query,reference}.rds
  first_level <- list.dirs(base_dir, recursive = FALSE, full.names = TRUE)
  first_level <- first_level[!grepl("\\.DS_Store", first_level)]

  rep_dirs <- first_level[grepl("/rep[0-9]+$", first_level)]

  if (length(rep_dirs) > 0) {
    dataset_dirs <- unlist(lapply(rep_dirs, function(rd) {
      list.dirs(rd, recursive = FALSE, full.names = TRUE)
    }))
  } else {
    dataset_dirs <- first_level
  }

  dataset_dirs <- dataset_dirs[!grepl("\\.DS_Store", dataset_dirs)]
  dataset_names <- basename(dataset_dirs)
  
  message_step("START", sprintf(
    "LABEL TRANSFER CLASSIFICATION PIPELINE\nFound %d datasets: %s",
    length(dataset_dirs),
    paste(dataset_names, collapse = ", ")
  ))
  
  # Process each dataset sequentially (can use parallel if needed later)
  message_step("RUN", "Processing datasets sequentially...")
  
  results <- list()
  for (i in seq_along(dataset_dirs)) {
    dataset_name <- dataset_names[i]
    dataset_dir <- dataset_dirs[i]
    
    message_step("DATASET", sprintf("[%d/%d] %s", i, length(dataset_dirs), dataset_name))
    
    result <- tryCatch(
      classify_dataset(dataset_dir, dataset_name, params, ncores = get_lt_n_cores(params)),
      error = function(e) {
        message_step("ERROR", sprintf("Dataset failed: %s", e$message))
        NULL
      }
    )
    
    results[[dataset_name]] <- result
  }
  
  # Summary
  n_success <- sum(!vapply(results, is.null, logical(1)))
  message_step("COMPLETE", sprintf(
    "Processed %d/%d datasets successfully",
    n_success, length(results)
  ))
  
  invisible(results)
}

# ============================================================================
# SCRIPT EXECUTION (when run from command line)
# ============================================================================

if (interactive() == FALSE) {
  # Script was called from command line
  args <- commandArgs(trailingOnly = TRUE)
  
  base_dir <- if (length(args) > 0) args[1] else NULL
  
  run_all_classifiers(base_dir = base_dir, params = params)
}
