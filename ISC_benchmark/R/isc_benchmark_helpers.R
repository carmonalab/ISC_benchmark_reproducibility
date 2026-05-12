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
source("R/assays_utils.R")
source("R/Metrics_benchmarking.R")

#' Get default gene blacklist (TCR, Immunoglobulins, Y-genes) from scTypeEval
get_default_blacklist <- function() {
  data("black_list", package = "scTypeEval", envir = environment())
  unlist(list(
    black_list$TCR,
    black_list$Immunoglobulins,
    black_list$Ygenes
  ))
}



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
    ident = ident_col
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
# BASELINE ISC (full dataset, unperturbed) CACHING
# ==========================================================================
#
# Design:
#   1. get_or_compute_baseline()   – computes/loads the full-dataset baseline
#      consistency dataframe (rate=1) once and caches it as baseline_isc_<ident>.rds.
#
#   2. baseline_for_task()         – relabels baseline df for tasks 1, 2, 5, 6.
#   3. baseline_for_Nct()          – relabels baseline df for task 3.
#   4. baseline_for_mergeCT()      – relabels baseline df for task 4.
#
#   5. Tasks 7/8: individual and merged subset scTypeEval objects are cached in
#      per-task output directories as subset_sc_* / merged_sc_* files.

#' Compute or load unified baseline consistency dataframe for all tasks
#'
#' Runs wrapper_scTypeEval once on full dataset and extracts consistency dataframe.
#' This single dataframe is used as baseline for all tasks (1-6).
#' Caches only the consistency dataframe, not the full scTypeEval object.
get_or_compute_baseline <- function(obj_prepared, config, cache_path) {
  if (file.exists(cache_path)) {
    message_step("BASELINE", sprintf("Loading cached baseline df from %s", basename(cache_path)))
    return(readRDS(cache_path))
  }

  message_step("BASELINE", "Computing full-dataset consistency baseline...")

  sc <- scTypeEval::create_scTypeEval(
    matrix = obj_prepared$count_matrix,
    metadata = obj_prepared$metadata,
    active_ident = obj_prepared$ident
  )

  sample_col  <- config$common$sample
  ndim        <- config$common$ndim
  reduction   <- config$common$reduction
  norm_method <- config$common$normalization_method
  min_samples <- config$common$min_samples
  min_cells   <- config$common$min_cells
  diss_method <- config$common$dissimilarity_method
  verbose_opt <- isTRUE(config$common$verbose)

  sc <- scTypeEval::wrapper_scTypeEval(
    sc,
    ident                = obj_prepared$ident,
    sample               = sample_col,
    gene_list            = NULL,
    reduction            = reduction,
    ndim                 = ndim,
    normalization_method = norm_method,
    dissimilarity_method = diss_method,
    min_samples          = min_samples,
    min_cells            = min_cells,
    verbose              = verbose_opt
  )

  baseline_df <- scTypeEval::get_consistency(sc, verbose = verbose_opt) |>
    dplyr::mutate(
      rate           = 1,
      rep            = 1,
      original_ident = obj_prepared$ident,
      task           = "Baseline"
    )

  saveRDS(baseline_df, cache_path)
  message_step("BASELINE", sprintf("Cached baseline df to %s", basename(cache_path)))
  baseline_df
}

#' Re-label a rate=1 baseline dataframe for tasks 1, 2, 5, 6
.baseline_task_labels <- list(
  missclassify = "Missclassification",
  SplitCelltype = "SplitCelltype",
  Nsamples     = "NSamples",
  NCell        = "NCell"
)

baseline_for_task <- function(baseline_df, task_name) {
  task_label <- .baseline_task_labels[[task_name]]
  if (is.null(task_label)) stop("Unknown baseline task: ", task_name)
  baseline_df$task <- task_label
  baseline_df$perturbed_ctype <- NA_character_
  baseline_df
}

#' Extract Nct baseline row from baseline dataframe (task 3)
#' Labels baseline dataframe row as task="Nct" with rate = all cell types joined by "-".
baseline_for_Nct <- function(baseline_df, all_cts) {
  all_cts_str <- paste(sort(as.character(all_cts)), collapse = "-")
  baseline_df |>
    dplyr::mutate(
      rate           = all_cts_str,
      rep            = NA_integer_,
      perturbed_ctype = NA_character_,
      task           = "Nct"
    )
}

#' Extract mergeCT baseline row from baseline dataframe (task 4)
#' Labels baseline dataframe row as task="mergeCT" with rate = original number of CTs.
baseline_for_mergeCT <- function(baseline_df, n_cts) {
  baseline_df |>
    dplyr::mutate(
      rate           = as.numeric(n_cts),
      rep            = NA_integer_,
      perturbed_ctype = NA_character_,
      task           = "mergeCT"
    )
}

#' Map batch/condition values to the dataset stem they originate from
#'
#' After merging datasets with add.cell.ids = dataset_stems, cell names carry a
#' "{stem}_" prefix. This function inverts that to give a batch → stem lookup.
create_stem_group_map <- function(metadata, dataset_stems, group_col) {
  if (is.null(dataset_stems) || length(dataset_stems) == 0 ||
      !group_col %in% colnames(metadata)) return(NULL)
  map <- list()
  cell_names <- rownames(metadata)
  for (stem in dataset_stems) {
    prefix      <- paste0(stem, "_")
    stem_cells  <- cell_names[startsWith(cell_names, prefix)]
    if (length(stem_cells) == 0) next
    groups <- unique(metadata[stem_cells, group_col, drop = TRUE])
    for (g in groups[!is.na(groups)]) {
      g_chr <- as.character(g)
      if (is.null(map[[g_chr]])) {
        map[[g_chr]] <- stem
      } else if (!identical(map[[g_chr]], stem)) {
        # Mark ambiguous group-to-stem mapping to avoid wrong baseline reuse.
        map[[g_chr]] <- NA_character_
      }
    }
  }
  map
}

#' Resolve the stems from a merged family that match a pair token set
#'
#' Used by tasks 7/8 to locate the two stems that belong to a valid pair.
#' Tokens are matched as case-insensitive substrings against dataset_stems.
resolve_pair_stems <- function(dataset_stems, tokens) {
  dataset_stems <- unique(trimws(unlist(strsplit(paste(dataset_stems, collapse = ","), ",", fixed = TRUE))))
  dataset_stems <- dataset_stems[nzchar(dataset_stems)]
  tokens <- unique(trimws(as.character(tokens)))
  tokens <- tokens[nzchar(tokens)]

  if (length(dataset_stems) == 0 || length(tokens) == 0) return(character(0))

  candidates <- dataset_stems
  for (token in tokens) {
    candidates <- candidates[grepl(token, candidates, fixed = TRUE, ignore.case = TRUE)]
  }
  unique(candidates)
}

#' Run scTypeEval on a cell subset with two-level caching (disk → memory)
#'
#' 1. Checks disk_cache_path first (persists across restarts; tasks 1-6 baseline
#'    files are valid candidates when dataset_stems are single files).
#' 2. Falls back to the in-memory cache_env for within-session deduplication.
#' 3. Computes fresh and saves to both levels if neither is available.
run_cached_subset_scTypeEval <- function(count_matrix,
                                         metadata,
                                         cell_idx,
                                         ident,
                                         config,
                                         cache_env,
                                         cache_key,
                                         cache_label,
                                         disk_cache_path = NULL) {
  # --- Level 1: in-memory cache ---
  if (exists(cache_key, envir = cache_env, inherits = FALSE)) {
    message(sprintf("    [mem cache] Reusing ISC for %s", cache_label))
    return(get(cache_key, envir = cache_env, inherits = FALSE))
  }

  # --- Level 2: disk cache for subset_sc_* / merged_sc_* objects ---
  if (!is.null(disk_cache_path) && file.exists(disk_cache_path)) {
    message(sprintf("    [disk cache] Loading ISC for %s from %s",
                    cache_label, basename(disk_cache_path)))
    result <- readRDS(disk_cache_path)
    assign(cache_key, result, envir = cache_env)
    return(result)
  }

  # --- Compute ---
  sample_col  <- config$common$sample
  ndim        <- config$common$ndim
  reduction   <- config$common$reduction
  norm_method <- config$common$normalization_method
  min_samples <- config$common$min_samples
  min_cells   <- config$common$min_cells
  diss_method <- config$common$dissimilarity_method
  verbose_opt <- isTRUE(config$common$verbose)

  sc <- scTypeEval::create_scTypeEval(
    matrix   = count_matrix[, cell_idx, drop = FALSE],
    metadata = metadata[cell_idx, , drop = FALSE],
    active_ident = ident
  )

  result <- scTypeEval::wrapper_scTypeEval(
    sc,
    ident                = ident,
    sample               = sample_col,
    gene_list            = NULL,
    reduction            = reduction,
    ndim                 = ndim,
    normalization_method = norm_method,
    dissimilarity_method = diss_method,
    min_samples          = min_samples,
    min_cells            = min_cells,
    verbose              = verbose_opt
  )

  # Save to disk if path given
  if (!is.null(disk_cache_path)) {
    saveRDS(result, disk_cache_path)
  }

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
  
  specs <- read.csv(specs_path,
                    stringsAsFactors = FALSE,
                      header = TRUE)
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

  batch_comp_col <- "Batch comparison"
  if (!batch_comp_col %in% names(specs)) {
    stop("Column not found in specs: ", batch_comp_col)
  }
  
  # Filter for batch comparison datasets
  batch_specs <- specs %>%
    filter(.data[[batch_comp_col]] == "yes") %>%
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

  pert_comp_col <- "Perturbation comparison"
  if (!pert_comp_col %in% names(specs)) {
    stop("Column not found in specs: ", pert_comp_col)
  }
  
  # Filter for perturbation comparison datasets
  pert_specs <- specs %>%
    filter(.data[[pert_comp_col]] == "yes") %>%
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
    # Disable internal wr_* file writing; orchestrator saves final baseline-inclusive result.
    list(dir = NULL)
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
    list(dir = NULL)
  )
  
  wr <- do.call(wr_split_cell_type, params)

  if (!is.null(baseline_df)) {
    wr <- rbind(baseline_for_task(baseline_df, "SplitCelltype"), wr)
  }

  wr
}

#' Run Task 3: Robustness to annotation granularity (coarse vs fine)
#' @return Result list from wr_nct()
run_task_Nct <- function(obj_prepared, config, task_config, output_dir, baseline_df = NULL) {
  message("Running Task 3: Robustness to annotation granularity")

  if (!is.null(baseline_df)) {
    task_config$run_baseline <- FALSE
    message("  [baseline reuse] Skipping full-CT run; will prepend cached baseline")
  }

  params <- c(
    obj_prepared,
    config$common,
    task_config,
    list(dir = NULL)
  )

  wr <- do.call(wr_nct, params)

  if (!is.null(baseline_df)) {
    all_cts <- unique(obj_prepared$metadata[[obj_prepared$ident]])
    all_cts <- all_cts[!is.na(all_cts)]
    bl <- baseline_for_Nct(baseline_df, all_cts)
    wr <- rbind(bl, wr)
  }

  wr
}

#' Run Task 4: Robustness to cellular complexity (high vs low variability)
run_task_cellular_complexity <- function(obj_prepared, config, task_config, output_dir, baseline_df = NULL) {
  message("Running Task 4: Robustness to cellular complexity")

  if (!is.null(baseline_df)) {
    task_config$run_original <- FALSE
    message("  [baseline reuse] Skipping original-annotation run; will prepend cached baseline")
  }

  params <- c(
    obj_prepared,
    config$common,
    task_config,
    list(dir = NULL)
  )

  wr <- do.call(wr_merge_ct, params)

  if (!is.null(baseline_df)) {
    n_cts <- length(unique(obj_prepared$metadata[[obj_prepared$ident]]))
    bl <- baseline_for_mergeCT(baseline_df, n_cts)
    wr <- rbind(bl, wr)
  }

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
    list(dir = NULL)
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
    list(dir = NULL)
  )
  
  wr <- do.call(wr_ncell, params)

  if (!is.null(baseline_df)) {
    wr <- rbind(baseline_for_task(baseline_df, "NCell"), wr)
  }

  wr
}

#' Run Task 7: Robustness to batch effects (systematic technical differences)
#'
#' For individual batches, this function first tries to reuse per-stem baseline
#' consistency files (`baseline_isc_<ident>.rds`) from tasks 1-6. If unavailable,
#' it computes and caches subset scTypeEval objects under task output directory.
#' Merged (pair) objects are always computed from the merged subset.
#'
#' @param results_root  Root results directory (config$output$dir).
#' @param dataset_stems Character vector of stems that were merged into the
#'   current merged object.  Used to map batch values back to stem names.
run_task_batch_effects <- function(obj_prepared, config, task_config, output_dir,
                                   specs_path = NULL,
                                   results_root = NULL,
                                   dataset_stems = NULL) {
  message("Running Task 7: Robustness to batch effects")

  metadata    <- obj_prepared$metadata
  count_matrix <- obj_prepared$count_matrix
  ident       <- obj_prepared$ident
  verbose_opt <- isTRUE(config$common$verbose)

  # Load specs to identify valid batch pairs
  if (is.null(specs_path)) {
    specs_path <- file.path(dirname(dirname(getwd())), "data_processing",
                            "config", "specs_datasets.csv")
  }

  specs       <- load_dataset_specs(specs_path)
  batch_pairs <- get_batch_pairs(specs, batch_col)

  if (is.null(batch_pairs) || length(batch_pairs) == 0) {
    message("  No valid batch pairs found in specs_datasets.csv")
    return(NULL)
  }

  message(sprintf("  Found %d batch pair(s) for comparison", length(batch_pairs)))

  dataset_stems <- unique(trimws(unlist(strsplit(paste(dataset_stems, collapse = ","), ",", fixed = TRUE))))
  dataset_stems <- dataset_stems[nzchar(dataset_stems)]

  tidy_results     <- data.frame()
  single_isc_cache <- new.env(parent = emptyenv())

  for (pair_idx in seq_along(batch_pairs)) {
    tryCatch({
      pair       <- batch_pairs[[pair_idx]]
      batch1     <- pair$batch1
      batch2     <- pair$batch2
      pair_name  <- pair$pair_name
      dataset_ref <- pair$dataset

      batch1_stems <- resolve_pair_stems(dataset_stems, c(batch1, pair$condition))
      batch2_stems <- resolve_pair_stems(dataset_stems, c(batch2, pair$condition))

      if (length(batch1_stems) != 1 || length(batch2_stems) != 1) {
        message(sprintf("    WARNING: Could not resolve unique stems for pair %s (%s -> %s; %s -> %s); skipping",
                        pair_name,
                        batch1, paste(batch1_stems, collapse = ","),
                        batch2, paste(batch2_stems, collapse = ",")))
      } else {
        batch1_stem <- batch1_stems[[1]]
        batch2_stem <- batch2_stems[[1]]

        message(sprintf("  Processing batch pair %d/%d: %s (%s) + %s (%s)",
                        pair_idx, length(batch_pairs), batch1, batch1_stem, batch2, batch2_stem))

        resolve_baseline_df_path <- function(stem) {
          if (is.null(stem) || is.na(stem) || is.null(results_root)) return(NULL)
          file.path(results_root, stem, paste0("baseline_isc_", ident, ".rds"))
        }

        # --- resolve disk cache paths for individual batches ---
        resolve_disk_path <- function(stem) {
          file.path(output_dir, paste0("subset_sc_", stem, "_", ident, ".rds"))
        }

        get_single_batch_consistency <- function(stem, batch_cells) {
          baseline_path <- resolve_baseline_df_path(stem)
          if (!is.null(baseline_path) && file.exists(baseline_path)) {
            message(sprintf("    [baseline reuse] Loading baseline consistency for batch '%s' from %s",
                            stem, baseline_path))
            return(readRDS(baseline_path))
          }

          isc_batch <- run_cached_subset_scTypeEval(
            count_matrix    = count_matrix,
            metadata        = metadata,
            cell_idx        = batch_cells,
            ident           = ident,
            config          = config,
            cache_env       = single_isc_cache,
            cache_key       = paste0("batch:", stem),
            cache_label     = sprintf("batch '%s'", stem),
            disk_cache_path = resolve_disk_path(stem)
          )
          scTypeEval::get_consistency(isc_batch, verbose = verbose_opt)
        }

        batch1_cells <- which(startsWith(rownames(metadata), paste0(batch1_stem, "_")))
        batch2_cells <- which(startsWith(rownames(metadata), paste0(batch2_stem, "_")))

        batch1_tidy <- get_single_batch_consistency(batch1_stem, batch1_cells) |>
          dplyr::mutate(batch = batch1, dataset = dataset_ref, ident = ident)

        batch2_tidy <- get_single_batch_consistency(batch2_stem, batch2_cells) |>
          dplyr::mutate(batch = batch2, dataset = dataset_ref, ident = ident)

        # --- merged object: always compute fresh ---
        combined_cells <- c(batch1_cells, batch2_cells)
        message(sprintf("    Merging pair: %s + %s", batch1_stem, batch2_stem))
        isc_combined <- run_cached_subset_scTypeEval(
          count_matrix    = count_matrix,
          metadata        = metadata,
          cell_idx        = combined_cells,
          ident           = ident,
          config          = config,
          cache_env       = single_isc_cache,
          cache_key       = paste0("merged:", pair_name),
          cache_label     = sprintf("merged '%s'", pair_name),
          disk_cache_path = file.path(output_dir,
                                      paste0("merged_sc_", pair_name, "_", ident, ".rds"))
        )

        combined_tidy <- scTypeEval::get_consistency(isc_combined, verbose = verbose_opt) |>
          dplyr::mutate(batch = pair_name, dataset = dataset_ref, ident = ident)

        pair_tidy  <- rbind(batch1_tidy, batch2_tidy, combined_tidy)
        pair_scores <- get_ratio(pair_name, pair_tidy, col = "batch")
        tidy_results <- rbind(tidy_results, pair_scores)
      }
    }, error = function(e) {
      message(sprintf("    ERROR: Batch pair %d failed (%s). Continuing with next pair.",
                      pair_idx, e$message))
      NULL
    })
  }

  if (nrow(tidy_results) == 0) {
    message("  WARNING: No valid batch comparisons could be computed")
    return(NULL)
  }

  message(sprintf("  Computed degradation scores for %d batch pair(s)", length(batch_pairs)))
  tidy_results
}

#' Run Task 8: Robustness to biological perturbations (condition-driven changes)
#'
#' For individual conditions, this function first tries to reuse per-stem baseline
#' consistency files (`baseline_isc_<ident>.rds`) from tasks 1-6. If unavailable,
#' it computes and caches subset scTypeEval objects under task output directory.
#' Merged (pair) objects are always computed from the merged subset.
#'
#' @param results_root  Root results directory (config$output$dir).
#' @param dataset_stems Character vector of stems merged into this object.
run_task_biological_perturbations <- function(obj_prepared, config, task_config, output_dir,
                                              specs_path = NULL,
                                              results_root = NULL,
                                              dataset_stems = NULL) {
  message("Running Task 8: Robustness to biological perturbations")

  metadata     <- obj_prepared$metadata
  count_matrix <- obj_prepared$count_matrix
  condition_col <- task_config$condition_col
  ident        <- obj_prepared$ident
  verbose_opt <- isTRUE(config$common$verbose)

  if (is.null(specs_path)) {
    specs_path <- file.path(dirname(dirname(getwd())), "data_processing",
                            "config", "specs_datasets.csv")
  }

  specs           <- load_dataset_specs(specs_path)
  condition_pairs <- get_perturbation_pairs(specs, batch_col = "Batch")

  if (is.null(condition_pairs) || length(condition_pairs) == 0) {
    message("  No valid condition pairs found in specs_datasets.csv")
    return(NULL)
  }

  message(sprintf("  Found %d condition pair(s) for comparison", length(condition_pairs)))

  dataset_stems <- unique(trimws(unlist(strsplit(paste(dataset_stems, collapse = ","), ",", fixed = TRUE))))
  dataset_stems <- dataset_stems[nzchar(dataset_stems)]

  tidy_results     <- data.frame()
  single_isc_cache <- new.env(parent = emptyenv())

  for (pair_idx in seq_along(condition_pairs)) {
    tryCatch({
      pair        <- condition_pairs[[pair_idx]]
      cond1       <- pair$condition1
      cond2       <- pair$condition2
      pair_name   <- pair$pair_name
      dataset_ref  <- pair$dataset

      cond1_stems <- resolve_pair_stems(dataset_stems, c(pair$batch, cond1))
      cond2_stems <- resolve_pair_stems(dataset_stems, c(pair$batch, cond2))

      if (length(cond1_stems) != 1 || length(cond2_stems) != 1) {
        message(sprintf("    WARNING: Could not resolve unique stems for pair %s (%s -> %s; %s -> %s); skipping",
                        pair_name,
                        cond1, paste(cond1_stems, collapse = ","),
                        cond2, paste(cond2_stems, collapse = ",")))
      } else {
        cond1_stem <- cond1_stems[[1]]
        cond2_stem <- cond2_stems[[1]]

        message(sprintf("  Processing condition pair %d/%d: %s (%s) + %s (%s)",
                        pair_idx, length(condition_pairs), cond1, cond1_stem, cond2, cond2_stem))

        resolve_baseline_df_path <- function(stem) {
          if (is.null(stem) || is.na(stem) || is.null(results_root)) return(NULL)
          file.path(results_root, stem, paste0("baseline_isc_", ident, ".rds"))
        }

        resolve_disk_path <- function(stem) {
          file.path(output_dir, paste0("subset_sc_", stem, "_", ident, ".rds"))
        }

        get_single_condition_consistency <- function(stem, cond_cells) {
          baseline_path <- resolve_baseline_df_path(stem)
          if (!is.null(baseline_path) && file.exists(baseline_path)) {
            message(sprintf("    [baseline reuse] Loading baseline consistency for condition '%s' from %s",
                            stem, baseline_path))
            return(readRDS(baseline_path))
          }

          isc_cond <- run_cached_subset_scTypeEval(
            count_matrix    = count_matrix,
            metadata        = metadata,
            cell_idx        = cond_cells,
            ident           = ident,
            config          = config,
            cache_env       = single_isc_cache,
            cache_key       = paste0("cond:", stem),
            cache_label     = sprintf("condition '%s'", stem),
            disk_cache_path = resolve_disk_path(stem)
          )
          scTypeEval::get_consistency(isc_cond, verbose = verbose_opt)
        }

        cond1_cells <- which(startsWith(rownames(metadata), paste0(cond1_stem, "_")))
        cond2_cells <- which(startsWith(rownames(metadata), paste0(cond2_stem, "_")))

        cond1_tidy <- get_single_condition_consistency(cond1_stem, cond1_cells) |>
          dplyr::mutate(condition = cond1, dataset = dataset_ref, ident = ident)

        cond2_tidy <- get_single_condition_consistency(cond2_stem, cond2_cells) |>
          dplyr::mutate(condition = cond2, dataset = dataset_ref, ident = ident)

        # --- merged object: always compute fresh ---
        combined_cells <- c(cond1_cells, cond2_cells)
        message(sprintf("    Merging pair: %s + %s", cond1_stem, cond2_stem))
        isc_combined <- run_cached_subset_scTypeEval(
          count_matrix    = count_matrix,
          metadata        = metadata,
          cell_idx        = combined_cells,
          ident           = ident,
          config          = config,
          cache_env       = single_isc_cache,
          cache_key       = paste0("merged:", pair_name),
          cache_label     = sprintf("merged '%s'", pair_name),
          disk_cache_path = file.path(output_dir,
                                      paste0("merged_sc_", pair_name, "_", ident, ".rds"))
        )

        combined_tidy <- scTypeEval::get_consistency(isc_combined, verbose = verbose_opt) |>
          dplyr::mutate(condition = pair_name, dataset = dataset_ref, ident = ident)

        pair_tidy   <- rbind(cond1_tidy, cond2_tidy, combined_tidy)
        pair_scores <- get_ratio(pair_name, pair_tidy, col = "condition")
        tidy_results <- rbind(tidy_results, pair_scores)
      }
    }, error = function(e) {
      message(sprintf("    ERROR: Condition pair %d failed (%s). Continuing with next pair.",
                      pair_idx, e$message))
      NULL
    })
  }

  if (nrow(tidy_results) == 0) {
    message("  WARNING: No valid condition comparisons could be computed")
    return(NULL)
  }

  message(sprintf("  Computed degradation scores for %d condition pair(s)", length(condition_pairs)))
  tidy_results
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
  
  # Save the per-task metrics file; this also serves as the resume marker.
  metrics_file <- file.path(
    output_dir,
    sprintf("%s_%s_%s_metrics.rds", dataset_id, task_name, ident)
  )
  saveRDS(results, metrics_file)
  message(sprintf("Saved metrics: %s", metrics_file))
  
  # Save full wr_* object if requested (and provided)
  save_wr_flag <- isTRUE(save_wr)
  if (save_wr_flag && !is.null(wr_object)) {
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
