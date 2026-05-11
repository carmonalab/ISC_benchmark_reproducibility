#!/usr/bin/env Rscript
# ISC_benchmark/R/isc_benchmark_helpers.R
#
# Helper functions for ISC benchmark reproducibility
# Wraps scTypeEval wr_* functions and provides metric computation
#
# Key functions:
# - prepare_scTypeEval_object() - Convert Seurat to scTypeEval format
# - run_isc_task() - Execute single task with error handling
# - extract_task_metrics() - Compute consistency metrics from wr_* results
# - save_task_results() - Persist results with metadata

suppressPackageStartupMessages({
  library(scTypeEval)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(yaml)
  library(BiocParallel)
})

# Source benchmarking wrappers (wr_* functions) from local vendored scripts.
source("R/01_assays_utils.R")
source("R/01_Metrics_benchmarking.R")



# ============================================================================
# DATA PREPARATION FOR scTypeEval
# ============================================================================

#' Prepare Seurat object for scTypeEval processing
#'
#' Converts Seurat object to scTypeEval's internal format with proper
#' normalization, PCA reduction, and metadata filtering
#'
#' @param obj Seurat object
#' @param ident_col Cell type annotation column
#' @param config Configuration list (from YAML)
# ============================================================================
# BASELINE ISC (rate = 1 / single, unperturbed) CACHING
# ============================================================================

#'
#' @return List with count_matrix, metadata, and metadata_filtered
prepare_scTypeEval_object <- function(obj, ident_col, config) {
  message(sprintf("[%s] Preparing object for scTypeEval", format(Sys.time(), "%H:%M:%S")))
  
  # Extract metadata and filter out NA identities
  metadata <- obj@meta.data
  if (!ident_col %in% colnames(metadata)) {
    stop("Cell type column '", ident_col, "' not found in metadata")
  }
  
  # Filter rows with NA cell type
  valid_cells <- !is.na(metadata[[ident_col]])
  metadata_filtered <- metadata[valid_cells, , drop = FALSE]
  count_matrix <- obj@assays$RNA$counts[, valid_cells]
  
  message(sprintf("  Cells: %d → %d (removed %d with NA %s)",
                  ncol(obj), ncol(count_matrix), sum(!valid_cells), ident_col))
  
  list(
    count_matrix = count_matrix,
    metadata = metadata_filtered,
    ident = ident_col,
    dir = NULL  # Will be set per task
  )
}

#' Load and merge multiple processed ISC datasets
#'
#' @param dataset_stems Character vector of processed dataset stems
#' @param config Configuration list
#'
#' @return Merged Seurat object containing all requested datasets
load_and_merge_processed_datasets <- function(dataset_stems, config) {
  dataset_stems <- unique(trimws(unlist(strsplit(paste(dataset_stems, collapse = ","), ",", fixed = TRUE))))
  dataset_stems <- dataset_stems[nzchar(dataset_stems)]

  if (length(dataset_stems) == 0) {
    stop("No dataset stems supplied for merging")
  }

  dataset_files <- file.path(config$processed_data_dir, paste0(dataset_stems, ".rds"))
  missing_files <- dataset_stems[!file.exists(dataset_files)]
  if (length(missing_files) > 0) {
    stop("Missing processed dataset(s): ", paste(missing_files, collapse = ", "))
  }

  objects <- lapply(dataset_files, readRDS)
  if (length(objects) == 1) {
    return(objects[[1]])
  }

  merged <- merge(
    x = objects[[1]],
    y = objects[-1],
    add.cell.ids = dataset_stems
  )
  merged
}

# ==========================================================================
# BASELINE ISC (rate = 1) CACHING
# ==========================================================================

.baseline_task_labels <- list(
  missclassify = "Missclassification",
  SplitCelltype = "SplitCelltype",
  Nsamples = "NSamples",
  NCell = "NCell"
)

#' Compute or load the unperturbed baseline ISC result
#'
#' This runs wr_missclasify once at rate = 1 and caches the output on disk so
#' tasks 1, 2, 5, and 6 can reuse the same baseline for a given dataset/ident.
get_or_compute_baseline <- function(obj_prepared, config, cache_path) {
  if (file.exists(cache_path)) {
    message_step("BASELINE", sprintf("Loading cached baseline from %s", basename(cache_path)))
    return(readRDS(cache_path))
  }

  message_step("BASELINE", "Computing unperturbed baseline (rate = 1)")

  params <- c(
    obj_prepared,
    config$common,
    list(
      rates = 1,
      replicates = 1,
      dir = NULL,
      save_plots = FALSE
    )
  )

  baseline_df <- do.call(wr_missclasify, params)
  saveRDS(baseline_df, cache_path)

  message_step("BASELINE", sprintf("Cached baseline to %s", basename(cache_path)))
  baseline_df
}

#' Re-label a baseline dataframe for a specific task
baseline_for_task <- function(baseline_df, task_name) {
  task_label <- .baseline_task_labels[[task_name]]
  if (is.null(task_label)) {
    stop("Unknown baseline task: ", task_name)
  }

  baseline_df$task <- task_label
  baseline_df
}

#' Run scTypeEval on a cell subset with in-memory caching
#'
#' The same single-batch or single-condition subset can appear in multiple
#' comparisons inside tasks 7 and 8. This helper makes those repeated calls
#' return the exact same scTypeEval result without recomputation.
run_cached_subset_scTypeEval <- function(count_matrix,
                                         metadata,
                                         cell_idx,
                                         ident,
                                         config,
                                         cache_env,
                                         cache_key,
                                         cache_label) {
  if (exists(cache_key, envir = cache_env, inherits = FALSE)) {
    message(sprintf("    [cache hit] Reusing ISC for %s", cache_label))
    return(get(cache_key, envir = cache_env, inherits = FALSE))
  }

  result <- scTypeEval::wrapper_scTypeEval(
    count_matrix = count_matrix[, cell_idx, drop = FALSE],
    metadata = metadata[cell_idx, , drop = FALSE],
    ident = ident,
    sample = config$common$sample_col,
    dissimilarity_method = config$common$dissimilarity_method,
    reduction = config$common$reduction,
    ndim = config$common$ndim,
    verbose = FALSE
  )

  assign(cache_key, result, envir = cache_env)
  result
}

# ============================================================================
# BATCH EFFECTS & PERTURBATIONS UTILITY FUNCTIONS
# ============================================================================

#' Load dataset specifications for batch and perturbation comparisons
#'
#' Reads specs_datasets.csv to identify which datasets can be used for
#' batch effects and perturbation comparison tasks.
#'
#' @param specs_path Path to specs_datasets.csv file
#'
#' @return Data frame with dataset specifications
load_dataset_specs <- function(specs_path) {
  if (!file.exists(specs_path)) {
    warning("Dataset specs file not found: ", specs_path)
    return(NULL)
  }
  
  specs <- read.csv(specs_path, stringsAsFactors = FALSE)
  specs
}

#' Create batch comparison pairs from dataset registry
#'
#' Identifies datasets that can be paired for batch effects comparison:
#' - Must have identical Dataset reference, Annotation, and Condition
#' - Must have different Batch values
#' - Must have "yes" in "Batch comparison" column
#'
#' @param obj_prepared Prepared scTypeEval object
#' @param specs Dataset specifications data frame
#' @param metadata Object metadata
#' @param batch_col Name of batch column in metadata
#'
#' @return List of batch pairs with metadata
#' @details
#' Example pairs from specs_datasets.csv:
#' - Cambridge Covid + Ncl Covid (same dataset, annotation, condition; different batch)
#' - CRC.SG1 Tumor + KUL3 Tumor (same dataset, annotation, condition; different batch)
get_batch_pairs <- function(specs, batch_col = "batch") {
  if (is.null(specs)) return(NULL)
  
  # Filter for batch comparison datasets
  batch_specs <- specs %>%
    filter(get(colnames(specs)[10]) == "yes") %>%  # "Batch comparison" column
    mutate(dataset_key = sprintf("%s_%s_%s",
                                 `Dataset reference`,
                                 Annotation,
                                 Condition))
  
  if (nrow(batch_specs) == 0) return(NULL)
  
  # Group by dataset_key to find pairs with same annotation+condition but different batches
  pairs <- batch_specs %>%
    group_by(dataset_key) %>%
    filter(n_distinct(Batch) > 1) %>%
    summarise(
      dataset_ref = first(`Dataset reference`),
      annotation = first(Annotation),
      condition = first(Condition),
      batches = list(unique(Batch)),
      .groups = "drop"
    )
  
  # Expand pairs to all combinations
  pair_list <- list()
  for (i in seq_len(nrow(pairs))) {
    batches <- pairs$batches[[i]]
    if (length(batches) >= 2) {
      combos <- combn(batches, 2, simplify = FALSE)
      for (combo in combos) {
        pair_list[[length(pair_list) + 1]] <- list(
          dataset = pairs$dataset_ref[i],
          annotation = pairs$annotation[i],
          condition = pairs$condition[i],
          batch1 = combo[1],
          batch2 = combo[2],
          pair_name = sprintf("%s-%s",  # Simple format for get_ratio()
                             combo[1],
                             combo[2])
        )
      }
    }
  }
  
  pair_list
}

#' Create perturbation comparison pairs from dataset registry
#'
#' Identifies datasets that can be paired for biological perturbation comparison:
#' - Must have identical Dataset reference, Annotation, and Batch
#' - Must have different Condition values
#' - Must have "yes" in "Perturbation comparison" column
#'
#' @param specs Dataset specifications data frame
#' @param batch_col Name of batch column in metadata
#'
#' @return List of condition pairs with metadata
#' @details
#' Example pairs from specs_datasets.csv:
#' - Cambridge Covid + Cambridge Healthy (same dataset, annotation, batch; different condition)
#' - CRC.SG1 Normal + CRC.SG1 Tumor (same dataset, annotation, batch; different condition)
get_perturbation_pairs <- function(specs, batch_col = "batch") {
  if (is.null(specs)) return(NULL)
  
  # Filter for perturbation comparison datasets
  pert_specs <- specs %>%
    filter(get(colnames(specs)[11]) == "yes") %>%  # "Perturbation comparison" column
    mutate(dataset_key = sprintf("%s_%s_%s",
                                 `Dataset reference`,
                                 Annotation,
                                 Batch))
  
  if (nrow(pert_specs) == 0) return(NULL)
  
  # Group by dataset_key to find pairs with same annotation+batch but different conditions
  pairs <- pert_specs %>%
    group_by(dataset_key) %>%
    filter(n_distinct(Condition) > 1, Condition != "") %>%  # Exclude empty conditions
    summarise(
      dataset_ref = first(`Dataset reference`),
      annotation = first(Annotation),
      batch = first(Batch),
      conditions = list(unique(Condition)),
      .groups = "drop"
    )
  
  if (nrow(pairs) == 0) return(NULL)
  
  # Expand pairs to all combinations
  pair_list <- list()
  for (i in seq_len(nrow(pairs))) {
    conds <- pairs$conditions[[i]]
    if (length(conds) >= 2) {
      combos <- combn(conds, 2, simplify = FALSE)
      for (combo in combos) {
        pair_list[[length(pair_list) + 1]] <- list(
          dataset = pairs$dataset_ref[i],
          annotation = pairs$annotation[i],
          batch = pairs$batch[i],
          condition1 = combo[1],
          condition2 = combo[2],
          pair_name = sprintf("%s-%s",  # Simple format for get_ratio()
                             combo[1],
                             combo[2])
        )
      }
    }
  }
  
  pair_list
}

#' Compute degradation score for batch/perturbation comparison
#'
#' Calculates how much ISC degrades when combining single batches/conditions
#' vs keeping them separate.
#'
#' Formula:
#' - drop = mean(single_measure) - combined_measure
#' - If drop < 0 (improvement): clip to 0
#' - score = 1 - drop (ranges 0-1, where 1=no degradation, 0=complete degradation)
#'
#' @param single_measures Vector of ISC scores for individual batches/conditions
#' @param combined_measure ISC score for combined batches/conditions
#'
#' @return Degradation score (0-1)
compute_degradation_score <- function(single_measures, combined_measure) {
  if (any(is.na(c(single_measures, combined_measure)))) {
    return(NA_real_)
  }
  
  single_mean <- mean(single_measures, na.rm = TRUE)
  drop <- single_mean - combined_measure
  
  # Clip improvements (negative drops) to 0 - only penalize degradations
  drop <- max(0, drop)
  
  # Score: 1 - drop→ensures 0 degradation gets score=1, complete degradation→score=0
  score <- 1 - drop
  score
}

#' Compute ratio/degradation score from multi-batch/condition
#'
#' Takes tidy results with separate rows for single and combined batches/conditions
#' and computes degradation scores using the formula from Consistency_metrics_benchmark.
#'
#' Input format required:
#' - dissimilarity_method: method used (e.g., "pearson", "euclidean")
#' - consistency.metric: consistency metric type (e.g., "ISC", "PCC")
#' - celltype: cell type identifier
#' - dataset: dataset ID
#' - ident: annotation/identity column name
#' - batch/condition: identifier (either single like "batch1" or multi like "batch1-batch2")
#' - measure: consistency value
#'
#' @param mb Multi-batch/condition identifier (e.g., "Cambridge-Ncl" for batch or "Covid-Healthy" for condition)
#' @param df Tidy data frame with measures (must have columns: dissimilarity_method, consistency.metric, celltype, dataset, ident, batch/condition, measure)
#' @param col Column name for grouping ("batch" or "condition")
#'
#' @return Data frame with degradation scores (method, metric, celltype, batch/condition, degradation_score)
get_ratio <- function(mb, df, col = "batch") {
  singles <- stringr::str_split(mb, "-", simplify = TRUE) %>% as.vector()
  
  group_cols <- c("dissimilarity_method", "consistency_metric", "celltype", "dataset", "ident")
  
  r <- df %>%
    filter(.data[[col]] %in% c(mb, singles)) %>%
    group_by(across(all_of(group_cols))) %>%
    summarise(
      multi_measure = measure[.data[[col]] == mb][1],
      single_mean = mean(measure[.data[[col]] %in% singles], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(!is.na(multi_measure), !is.na(single_mean)) %>%
    mutate(
      !!col := mb,
      drop = single_mean - multi_measure,
      # Set improvements (negatives) to 0
      score_celltype = if_else(drop < 0, 0, drop),
      degradation_score = 1 - score_celltype
    ) %>%
    select(all_of(group_cols), !!col, degradation_score)
  
  return(r)
}

# ============================================================================
# TASK EXECUTION WRAPPERS (call scTypeEval wr_* functions)
# ============================================================================

#' Run Task 1: Sensitivity to cell type signal degradation (label noise)
run_task_missclassify <- function(obj_prepared, config, task_config, output_dir,
                                   baseline_df = NULL) {
  message("Running Task 1: Sensitivity to cell type signal degradation (Label Noise)")
  
  if (!is.null(baseline_df)) {
    task_config$rates <- task_config$rates[task_config$rates != 1]
    message("  [baseline reuse] Skipping rate=1; will prepend cached baseline")
  }

  params <- c(
    obj_prepared,
    config$common,
    task_config,
    list(dir = output_dir)
  )
  
  wr <- do.call(wr_missclasify, params)

  if (!is.null(baseline_df)) {
    wr <- rbind(baseline_for_task(baseline_df, "missclassify"), wr)
  }

  wr
}

#' Run Task 2: Sensitivity to cell type over-partitioning (artificial subtypes)
run_task_SplitCelltype <- function(obj_prepared, config, task_config, output_dir,
                                    baseline_df = NULL) {
  message("Running Task 2: Sensitivity to cell type over-partitioning")
  
  if (!is.null(baseline_df)) {
    task_config$rates <- task_config$rates[task_config$rates != 1]
    message("  [baseline reuse] Skipping rate=1; will prepend cached baseline")
  }

  params <- c(
    obj_prepared,
    config$common,
    task_config,
    list(dir = output_dir)
  )
  
  wr <- do.call(wr_split_cell_type, params)

  if (!is.null(baseline_df)) {
    wr <- rbind(baseline_for_task(baseline_df, "SplitCelltype"), wr)
  }

  wr
}

#' Run Task 3: Robustness to annotation granularity (coarse vs fine)
#' @return Result list from scTypeEval:::wr.Nct()
run_task_Nct <- function(obj_prepared, config, task_config, output_dir) {
  message("Running Task 3: Robustness to annotation granularity")
  
  params <- c(
    obj_prepared,
    config$common,
    task_config,
    list(dir = output_dir)
  )
  
  wr <- do.call(wr_nct, params)
  wr
}

#' Run Task 4: Robustness to cellular complexity (high vs low variability)
run_task_cellular_complexity <- function(obj_prepared, config, task_config, output_dir) {
  message("Running Task 4: Robustness to cellular complexity")
 
   params <- c(
    obj_prepared,
    config$common,
    task_config,
    list(dir = output_dir)
  )
  
  wr <- do.call(wr_merge_ct, params)
  wr
}

#' Run Task 5: Robustness to dataset size - samples (varying sample count)
run_task_Nsamples <- function(obj_prepared, config, task_config, output_dir,
                               baseline_df = NULL) {
  message("Running Task 5: Robustness to dataset size (samples)")
  
  if (!is.null(baseline_df)) {
    task_config$rates <- task_config$rates[task_config$rates != 1]
    message("  [baseline reuse] Skipping rate=1; will prepend cached baseline")
  }

  params <- c(
    obj_prepared,
    config$common,
    task_config,
    list(dir = output_dir)
  )
  
  wr <- do.call(wr_nsamples, params)

  if (!is.null(baseline_df)) {
    wr <- rbind(baseline_for_task(baseline_df, "Nsamples"), wr)
  }

  wr
}

#' Run Task 6: Robustness to dataset size - cells (varying cells per cell type)
run_task_NCell <- function(obj_prepared, config, task_config, output_dir,
                            baseline_df = NULL) {
  message("Running Task 6: Robustness to dataset size (cells)")
  
  if (!is.null(baseline_df)) {
    task_config$rates <- task_config$rates[task_config$rates != 1]
    message("  [baseline reuse] Skipping rate=1; will prepend cached baseline")
  }

  params <- c(
    obj_prepared,
    config$common,
    task_config,
    list(dir = output_dir)
  )
  
  wr <- do.call(wr_ncell, params)

  if (!is.null(baseline_df)) {
    wr <- rbind(baseline_for_task(baseline_df, "NCell"), wr)
  }

  wr
}

#' Run Task 7: Robustness to batch effects (systematic technical differences)
#'
#' Compares ISC degradation when combining multiple batches vs keeping them separate.
#' Uses specs_datasets.csv to identify valid batch pairs:
#' - Same dataset reference + annotation + condition
#' - Different batches
#'
#' Returns tidy format results with batch column (single or multi-batch identifier)
#' and computes degradation scores using get_ratio()
#'
#' @param obj_prepared Prepared object from prepare_scTypeEval_object()
#' @param config Common configuration list
#' @param task_config Task-specific configuration from isc_benchmark_parameters.yaml
#' @param output_dir Output directory for this task
#' @param specs_path Path to specs_datasets.csv file (optional)
#'
#' @return Data frame with batch effect scores (method, metric, celltype, batch, degradation_score)
run_task_batch_effects <- function(obj_prepared, config, task_config, output_dir, 
                                   specs_path = NULL) {
  message("Running Task 7: Robustness to batch effects")
  
  metadata <- obj_prepared$metadata
  count_matrix <- obj_prepared$count_matrix
  batch_col <- task_config$batch_col
  
  if (!batch_col %in% colnames(metadata)) {
    message(sprintf("  WARNING: Batch column '%s' not found; returning NULL", batch_col))
    return(NULL)
  }
  
  # Load specs to identify valid batch pairs
  if (is.null(specs_path)) {
    specs_path <- file.path(dirname(dirname(getwd())), "data_processing", 
                            "config", "specs_datasets.csv")
  }
  
  specs <- load_dataset_specs(specs_path)
  batch_pairs <- get_batch_pairs(specs, batch_col)
  
  if (is.null(batch_pairs) || length(batch_pairs) == 0) {
    message("  No valid batch pairs found in specs_datasets.csv")
    return(NULL)
  }
  
  message(sprintf("  Found %d batch pair(s) for comparison", length(batch_pairs)))
  
  tidy_results <- data.frame()
  
  # Cache single-batch ISC results; the same batch can appear in multiple pairs.
  single_isc_cache <- new.env(parent = emptyenv())

  # Process each batch pair
  for (pair_idx in seq_along(batch_pairs)) {
    tryCatch({
      pair <- batch_pairs[[pair_idx]]
      batch1 <- pair$batch1
      batch2 <- pair$batch2
      pair_name <- pair$pair_name
      dataset_ref <- pair$dataset
      
      message(sprintf("  Processing batch pair %d/%d: %s + %s",
                      pair_idx, length(batch_pairs), batch1, batch2))
      
      # Get cells for each batch
      batch1_cells <- which(metadata[[batch_col]] == batch1)
      batch2_cells <- which(metadata[[batch_col]] == batch2)
      
      if (length(batch1_cells) < 50 || length(batch2_cells) < 50) {
        message(sprintf("    WARNING: Insufficient cells (%d, %d); skipping pair",
                        length(batch1_cells), length(batch2_cells)))
      } else {
        # Compute ISC for individual batches using the shared subset cache.
        isc_batch1 <- run_cached_subset_scTypeEval(
          count_matrix = count_matrix,
          metadata = metadata,
          cell_idx = batch1_cells,
          ident = obj_prepared$ident,
          config = config,
          cache_env = single_isc_cache,
          cache_key = paste0("batch:", batch1),
          cache_label = sprintf("batch '%s'", batch1)
        )

        isc_batch2 <- run_cached_subset_scTypeEval(
          count_matrix = count_matrix,
          metadata = metadata,
          cell_idx = batch2_cells,
          ident = obj_prepared$ident,
          config = config,
          cache_env = single_isc_cache,
          cache_key = paste0("batch:", batch2),
          cache_label = sprintf("batch '%s'", batch2)
        )
        
        # Compute ISC for combined batches
        combined_cells <- c(batch1_cells, batch2_cells)
        isc_combined <- scTypeEval::wrapper_scTypeEval(
          count_matrix = count_matrix[, combined_cells],
          metadata = metadata[combined_cells, , drop = FALSE],
          ident = obj_prepared$ident,
          sample = config$common$sample_col,
          dissimilarity_method = config$common$dissimilarity_method,
          reduction = config$common$reduction,
          ndim = config$common$ndim,
          verbose = FALSE
        )
        
        # Convert to tidy format for get_ratio()
        batch1_tidy <- scTypeEval::get_consistency(isc_batch1) %>%
          mutate(batch = batch1, dataset = dataset_ref, ident = obj_prepared$ident)
        
        batch2_tidy <- scTypeEval::get_consistency(isc_batch2) %>%
          mutate(batch = batch2, dataset = dataset_ref, ident = obj_prepared$ident)
        
        combined_tidy <- scTypeEval::get_consistency(isc_combined) %>%
          mutate(batch = pair_name, dataset = dataset_ref, ident = obj_prepared$ident)
        
        # Combine into single tidy frame for this pair
        pair_tidy <- rbind(batch1_tidy, batch2_tidy, combined_tidy)
        
        # Compute degradation score using get_ratio()
        pair_scores <- get_ratio(pair_name, pair_tidy, col = "batch")
        
        tidy_results <- rbind(tidy_results, pair_scores)
      }
    }, error = function(e) {
      message(sprintf("    ERROR: Batch pair %d failed (%s). Continuing with next pair.",
                      pair_idx, e$message))
      NULL
    })
  }
  
  # Combine all results
  if (nrow(tidy_results) == 0) {
    message("  WARNING: No valid batch comparisons could be computed")
    return(NULL)
  }
  
  message(sprintf("  ✓ Computed degradation scores for %d batch pair(s)", length(batch_pairs)))
  tidy_results
}

#' Run Task 8: Robustness to biological perturbations (condition-driven changes)
#'
#' Compares ISC degradation when combining multiple biological conditions vs keeping them separate.
#' Uses specs_datasets.csv to identify valid condition pairs:
#' - Same dataset reference + annotation + batch
#' - Different conditions
#'
#' Computes:
#' 1. ISC for each individual condition (baseline)
#' 2. ISC for combined conditions
#' 3. Degradation score: 1 - (mean(single_ISC) - combined_ISC)
#'
#' @param obj_prepared Prepared object from prepare_scTypeEval_object()
#' @param config Common configuration list
#' @param task_config Task-specific configuration from isc_benchmark_parameters.yaml
#' @param output_dir Output directory for this task
#' @param specs_path Path to specs_datasets.csv file (optional)
#'
#' @return Data frame with perturbation effect scores (method, metric, celltype, condition_pair, score)
run_task_biological_perturbations <- function(obj_prepared, config, task_config, output_dir,
                                              specs_path = NULL) {
  message("Running Task 8: Robustness to biological perturbations")
  
  metadata <- obj_prepared$metadata
  count_matrix <- obj_prepared$count_matrix
  condition_col <- task_config$condition_col
  
  if (!condition_col %in% colnames(metadata)) {
    message(sprintf("  WARNING: Condition column '%s' not found; returning NULL", condition_col))
    return(NULL)
  }
  
  # Load specs to identify valid condition pairs
  if (is.null(specs_path)) {
    specs_path <- file.path(dirname(dirname(getwd())), "data_processing", 
                            "config", "specs_datasets.csv")
  }
  
  specs <- load_dataset_specs(specs_path)
  condition_pairs <- get_perturbation_pairs(specs, batch_col = "Batch")
  
  if (is.null(condition_pairs) || length(condition_pairs) == 0) {
    message("  No valid condition pairs found in specs_datasets.csv")
    return(NULL)
  }
  
  message(sprintf("  Found %d condition pair(s) for comparison", length(condition_pairs)))
  
  tidy_results <- data.frame()
  
  # Cache single-condition ISC results; the same condition can appear in multiple pairs.
  single_isc_cache <- new.env(parent = emptyenv())

  # Process each condition pair
  for (pair_idx in seq_along(condition_pairs)) {
    tryCatch({
      pair <- condition_pairs[[pair_idx]]
      cond1 <- pair$condition1
      cond2 <- pair$condition2
      pair_name <- pair$pair_name
      dataset_ref <- pair$dataset
      
      message(sprintf("  Processing condition pair %d/%d: %s + %s",
                      pair_idx, length(condition_pairs), cond1, cond2))
      
      # Get cells for each condition
      cond1_cells <- which(metadata[[condition_col]] == cond1)
      cond2_cells <- which(metadata[[condition_col]] == cond2)
      
      if (length(cond1_cells) < 50 || length(cond2_cells) < 50) {
        message(sprintf("    WARNING: Insufficient cells (%d, %d); skipping pair",
                        length(cond1_cells), length(cond2_cells)))
      } else {
        # Compute ISC for individual conditions using the shared subset cache.
        isc_cond1 <- run_cached_subset_scTypeEval(
          count_matrix = count_matrix,
          metadata = metadata,
          cell_idx = cond1_cells,
          ident = obj_prepared$ident,
          config = config,
          cache_env = single_isc_cache,
          cache_key = paste0("cond:", cond1),
          cache_label = sprintf("condition '%s'", cond1)
        )

        isc_cond2 <- run_cached_subset_scTypeEval(
          count_matrix = count_matrix,
          metadata = metadata,
          cell_idx = cond2_cells,
          ident = obj_prepared$ident,
          config = config,
          cache_env = single_isc_cache,
          cache_key = paste0("cond:", cond2),
          cache_label = sprintf("condition '%s'", cond2)
        )
        
        # Compute ISC for combined conditions
        combined_cells <- c(cond1_cells, cond2_cells)
        isc_combined <- scTypeEval::wrapper_scTypeEval(
          count_matrix = count_matrix[, combined_cells],
          metadata = metadata[combined_cells, , drop = FALSE],
          ident = obj_prepared$ident,
          sample = config$common$sample_col,
          dissimilarity_method = config$common$dissimilarity_method,
          reduction = config$common$reduction,
          ndim = config$common$ndim,
          verbose = FALSE
        )
        
        # Convert to tidy format for get_ratio()
        cond1_tidy <- scTypeEval::get_consistency(isc_cond1) %>%
          mutate(condition = cond1, dataset = dataset_ref, ident = obj_prepared$ident)
        
        cond2_tidy <- scTypeEval::get_consistency(isc_cond2) %>%
          mutate(condition = cond2, dataset = dataset_ref, ident = obj_prepared$ident)
        
        combined_tidy <- scTypeEval::get_consistency(isc_combined) %>%
          mutate(condition = pair_name, dataset = dataset_ref, ident = obj_prepared$ident)
        
        # Combine into single tidy frame for this pair
        # Note: pair_name is like "Healthy-Tumor" from get_perturbation_pairs
        pair_tidy <- rbind(cond1_tidy, cond2_tidy, combined_tidy)
        
        # Compute degradation score using get_ratio()
        pair_scores <- get_ratio(pair_name, pair_tidy, col = "condition")
        
        tidy_results <- rbind(tidy_results, pair_scores)
      }
    }, error = function(e) {
      message(sprintf("    ERROR: Condition pair %d failed (%s). Continuing with next pair.",
                      pair_idx, e$message))
      NULL
    })
  }
  
  # Combine all results
  if (nrow(tidy_results) == 0) {
    message("  WARNING: No valid condition comparisons could be computed")
    return(NULL)
  }
  
  message(sprintf("  ✓ Computed degradation scores for %d condition pair(s)", length(condition_pairs)))
  tidy_results
}

# ============================================================================
# METRIC EXTRACTION FROM TASK RESULTS
# ============================================================================

#' Extract metrics from wr_* task results
#'
#' Converts scTypeEval wr_* output into tidy data frame for analysis
#'
#' @param wr_result Result from scTypeEval wr_* function
#' @param task_name Name of task (for labeling)
#' @param metric_type Type of expected degradation ("monotonic" or "constant")
#'
#' @return Data frame with metrics, rates, and degradation scores
extract_task_metrics <- function(wr_result, task_name, metric_type) {
  
  # Use scTypeEval's built-in plotting function to extract metrics
  # This returns a tidy data frame suitable for analysis
  metrics_df <- wr_assay_plot(
    wr_result,
    type = metric_type,
    return_df = TRUE
  )
  
  # Add task metadata
  metrics_df <- metrics_df %>%
    mutate(task = task_name, .before = 1) %>%
    mutate(metric_type = metric_type, .after = "task")
  
  metrics_df
}

# ============================================================================
# RESULT PERSISTENCE AND LOGGING
# ============================================================================

#' Save task results with metadata
#'
#' Persists task results, metrics, and session information
#'
#' @param results Task results (metrics data frame)
#' @param wr_object Full wr_* object (optional, large)
#' @param task_name Name of task
#' @param dataset_id Dataset identifier
#' @param ident Cell type annotation name
#' @param output_dir Output directory
#' @param config Configuration list (for snapshot)
#' @param save_wr Whether to save full wr_* object
#'
#' @return Path to saved results file
save_task_results <- function(results,
                              wr_object = NULL,
                              task_name,
                              dataset_id,
                              ident,
                              output_dir,
                              config,
                              save_wr = FALSE) {
  
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Save metrics
  metrics_file <- file.path(
    output_dir,
    sprintf("%s_%s_%s_metrics.rds", dataset_id, task_name, ident)
  )
  saveRDS(results, metrics_file)
  message(sprintf("Saved metrics: %s", metrics_file))
  
  # Save full wr_* object if requested (and provided)
  if (save_wr && !is.null(wr_object)) {
    wr_file <- file.path(
      output_dir,
      sprintf("%s_%s_%s_wrobj.rds", dataset_id, task_name, ident)
    )
    saveRDS(wr_object, wr_file)
    message(sprintf("Saved wr_* object: %s", wr_file))
  }
  
  # Save metadata with config snapshot
  metadata <- list(
    task = task_name,
    dataset = dataset_id,
    ident = ident,
    timestamp = Sys.time(),
    seed = config$seed,
    n_cores = config$n_cores,
    task_config = config[[paste0("task_", task_name)]],
    common_config = config$common,
    session_info = utils::sessionInfo()
  )
  
  config_file <- file.path(
    output_dir,
    sprintf("%s_%s_%s_metadata.yaml", dataset_id, task_name, ident)
  )
  write_yaml(metadata, config_file)
  message(sprintf("Saved metadata: %s", config_file))
  
  invisible(metrics_file)
}
