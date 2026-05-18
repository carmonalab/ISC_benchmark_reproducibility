# label_transfer_task/R/01_classifiers.R --- Classifier execution integration

# This module integrates the main classifier implementations from 01_run_classifiers.R
# and provides helper functions for targets workflow

source("R/00_utils.R")
source("R/01_run_classifiers.R")

# ============================================================================
# TARGETS-FRIENDLY WRAPPER: Run a single classifier on a dataset
# ============================================================================

#' Run single classifier on dataset (for targets workflow)
#'
#' @param dataset_id Dataset identifier
#' @param classifier_name Classifier name (must match function in classifiers.R)
#' @param data_dir Directory with query.rds and reference.rds
#' @param output_dir Directory to save results
#' @param seed Random seed for reproducibility
#'
#' @return Data frame with predictions and metadata
run_label_transfer_classifier_targets <- function(
    dataset_id,
    classifier_name,
    rep = 1,
    data_dir = NULL,
    output_dir = NULL,
    seed = NULL,
  ncores = 1,
  reference_dataset_id = NULL,
  query_dataset_id = NULL) {

  if (is.null(seed)) {
    stop("Missing required argument: seed (should be provided from config)")
  }

  set.seed(as.integer(seed))
  
  if (is.null(data_dir)) {
    data_dir <- proj_path("data/processed/label_transfer")
  }
  
  if (is.null(output_dir)) {
    output_dir <- lt_raw_results_dir()
  }
  
  dataset_dir <- file.path(data_dir, dataset_id)
  
  # Load data
  query_file <- file.path(dataset_dir, "query.rds")
  reference_file <- file.path(dataset_dir, "reference.rds")
  
  if (!file.exists(query_file) || !file.exists(reference_file)) {
    return(NULL)
  }
  
  query <- readRDS(query_file)
  reference <- readRDS(reference_file)
  
  # Load classifier functions
  local_env <- new.env(parent = environment())
  local_env$ncores <- ncores
  source(proj_path("label_transfer_task/classifiers/classifiers.R"), local = local_env)
  
  # Get classifier function
  clf_func_name <- paste0("classify_", classifier_name)
  if (!exists(clf_func_name, envir = local_env)) {
    warning("Classifier function not found: ", clf_func_name)
    return(NULL)
  }
  
  clf_func <- get(clf_func_name, envir = local_env)
  
  # Run classifier
  predictions <- tryCatch({
    clf_func(
      reference$counts,
      reference$metadata$cell_type,
      query$counts
    )
  }, error = function(e) {
    message(sprintf("[SKIP] Classifier '%s' failed on '%s' (rep %d): %s",
                    classifier_name, dataset_id, rep, conditionMessage(e)))
    NA
  })

  # Treat scalar NA or length mismatch as classifier failure and do not emit an output file.
  n_query <- ncol(query$counts)
  failed <- (length(predictions) == 1L && is.na(predictions)) || (length(predictions) != n_query)
  if (failed) {
    stop(sprintf(
      "Classifier '%s' failed on '%s' (rep %d): predictions invalid (length=%d, expected=%d)",
      classifier_name, dataset_id, rep, length(predictions), n_query
    ))
  }

  # Ensure predictions can be aligned downstream
  if (!all(is.na(predictions)) && length(predictions) == n_query) {
    names(predictions) <- colnames(query$counts)
  }
  
  # Create result data frame
  result_df <- query$metadata %>%
    as.data.frame() %>%
    tibble::rownames_to_column("cell_id") %>%
    mutate(
      dataset_id = dataset_id,
      classifier = classifier_name,
      replicate = rep,
      prediction = unname(predictions),
      accuracy = if (all(is.na(predictions))) NA else 
                 mean(predictions == cell_type, na.rm = TRUE),
      seed = seed
    )

  if (!is.null(reference_dataset_id)) {
    result_df$reference_dataset_id <- reference_dataset_id
  }
  if (!is.null(query_dataset_id)) {
    result_df$query_dataset_id <- query_dataset_id
  }
  
  # Save result
  output_file <- file.path(
    output_dir,
    sprintf("%s_%s_rep%d.rds", dataset_id, classifier_name, rep)
  )
  saveRDS(result_df, output_file)

  return(output_file)
}

# Keep legacy function name for compatibility
run_label_transfer_classifier <- run_label_transfer_classifier_targets

# ============================================================================
# HELPER: Get all classifiers list
# ============================================================================

get_available_classifiers <- function() {
  c(
    "SingleR",
    "LogisticRegression",
    "XGBoost",
    "MLP",
    "RandomForest",
    "SVM",
    "kNN",
    "NaiveBayes",
    "LDA",
    "SeuratTransfer",
    "DecisionTree",
    "scPred",
    "Random",
    "Ensemble"
  )
}

# ============================================================================
# ENSEMBLE: Majority voting after all classifiers complete
# ============================================================================

#' Run Ensemble meta-classifier (majority voting)
#'
#' Aggregates predictions from all individual classifiers (except Random baseline)
#' and uses majority voting to assign final labels.
#'
#' @param dataset_id Dataset identifier
#' @param rep Replicate number
#' @param results_dir Directory containing individual classifier results
#' @param output_dir Directory to save Ensemble results
#' @param seed Random seed
#'
#' @return Path to saved Ensemble result file
run_ensemble_classifier_targets <- function(
    dataset_id,
    rep = 1,
    results_dir = NULL,
    output_dir = NULL,
    seed = NULL) {

  if (is.null(results_dir)) results_dir <- lt_raw_results_dir()
  if (is.null(output_dir)) output_dir <- lt_raw_results_dir()
  if (is.null(seed)) seed <- 22
  
  set.seed(as.integer(seed))
  
  # Load classifier functions
  local_env <- new.env(parent = environment())
  source(proj_path("label_transfer_task/classifiers/classifiers.R"), local = local_env)
  
  # Find all classifier results for this dataset/rep, excluding Random and Ensemble
  pattern <- sprintf("%s_.*_rep%d\\.rds$", dataset_id, rep)
  all_files <- list.files(results_dir, pattern = pattern, full.names = TRUE)
  
  # Load all predictions
  predictions_list <- list()
  first_result <- NULL
  for (file in all_files) {
    result <- tryCatch(readRDS(file), error = function(e) NULL)
    if (!is.null(result) && nrow(result) > 0) {
      classifier_name <- basename(file) %>%
        stringr::str_replace(sprintf("%s_", dataset_id), "") %>%
        stringr::str_replace(sprintf("_rep%d\\.rds", rep), "")
      
      # Exclude Random baseline and previous Ensemble attempts
      if (!classifier_name %in% c("Random", "Ensemble")) {
        predictions_list[[classifier_name]] <- result$prediction
        if (is.null(first_result)) first_result <- result
      }
    }
  }
  
  if (length(predictions_list) == 0) {
    message("  [Ensemble] No individual classifier results found, skipping.")
    return(NA_character_)
  }
  
  # Call existing classify_Ensemble function
  message(sprintf("  [Ensemble] Aggregating %d classifiers...", length(predictions_list)))
  
  ensemble_pred <- tryCatch(
    get("classify_Ensemble", envir = local_env)(predictions_list),
    error = function(e) {
      message("  [Ensemble] Voting failed: ", e$message)
      NA
    }
  )
  
  # Create result dataframe matching first_result structure
  result_df <- first_result %>%
    dplyr::select(-all_of(c("classifier", "prediction", "accuracy", "seed"))) %>%
    dplyr::mutate(
      classifier = "Ensemble",
      prediction = unname(ensemble_pred),
      accuracy = if (all(is.na(ensemble_pred))) NA else 
                 mean(ensemble_pred == cell_type, na.rm = TRUE),
      seed = as.integer(seed)
    )
  
  # Save result
  output_file <- file.path(
    output_dir,
    sprintf("%s_Ensemble_rep%d.rds", dataset_id, rep)
  )
  saveRDS(result_df, output_file)
  message(sprintf("    ✓ Ensemble result saved to %s", basename(output_file)))
  
  return(output_file)
}

# ============================================================================
# ENSEMBLE: Majority voting for between-dataset transfers
# ============================================================================

#' Run Ensemble meta-classifier for between-dataset transfers
#'
#' Aggregates predictions from all classifiers (except Random and Ensemble) 
#' for a specific between-dataset transfer pair using the existing classify_Ensemble function.
#'
#' @param pair_id Between-dataset pair identifier (e.g., "source_target")
#' @param reference_dataset_id Reference dataset ID
#' @param query_dataset_id Query dataset ID
#' @param rep Replicate number
#' @param output_dir Directory to save Ensemble results
#'
#' @return Path to saved Ensemble result file
run_ensemble_classifier_between_datasets <- function(
    pair_id,
    reference_dataset_id,
    query_dataset_id,
    rep = 1,
    output_dir = NULL) {

  if (is.null(output_dir)) {
    output_dir <- lt_between_raw_results_dir()
  }

  # Load classifier functions
  local_env <- new.env(parent = environment())
  source(proj_path("label_transfer_task/classifiers/classifiers.R"), local = local_env)

  # Find all classifier results for this pair/rep, excluding Random and Ensemble
  pattern <- sprintf("%s_.*_rep%d\\.rds$", pair_id, rep)
  all_files <- list.files(output_dir, pattern = pattern, full.names = TRUE)
  
  # Load all predictions
  predictions_list <- list()
  first_result <- NULL
  
  for (file in all_files) {
    result <- tryCatch(readRDS(file), error = function(e) NULL)
    if (!is.null(result) && nrow(result) > 0) {
      classifier_name <- basename(file) %>%
        stringr::str_replace(sprintf("%s_", pair_id), "") %>%
        stringr::str_replace(sprintf("_rep%d\\.rds", rep), "")
      
      # Exclude Random baseline and previous Ensemble attempts
      if (!classifier_name %in% c("Random", "Ensemble")) {
        predictions_list[[classifier_name]] <- result$prediction
        if (is.null(first_result)) first_result <- result
      }
    }
  }

  if (length(predictions_list) == 0) {
    message("  [Ensemble] No individual classifier results found for pair ", pair_id, ", skipping.")
    return(NA_character_)
  }

  # Call existing classify_Ensemble function
  message(sprintf("  [Ensemble] Aggregating %d classifiers for pair %s...", 
                  length(predictions_list), pair_id))
  
  ensemble_pred <- tryCatch(
    get("classify_Ensemble", envir = local_env)(predictions_list),
    error = function(e) {
      message("  [Ensemble] Voting failed for pair ", pair_id, ": ", e$message)
      NA
    }
  )

  # Create result data frame (same structure as individual classifiers)
  result_df <- first_result %>%
    dplyr::select(-all_of(c("classifier", "prediction", "accuracy", "seed"))) %>%
    dplyr::mutate(
      classifier = "Ensemble",
      prediction = unname(ensemble_pred),
      accuracy = if (all(is.na(ensemble_pred))) NA else 
                 mean(ensemble_pred == cell_type, na.rm = TRUE),
      seed = NA_integer_,
      reference_dataset_id = reference_dataset_id,
      query_dataset_id = query_dataset_id
    )

  # Save result
  output_file <- file.path(
    output_dir,
    sprintf("%s_Ensemble_rep%d.rds", pair_id, rep)
  )
  saveRDS(result_df, output_file)
  message(sprintf("    ✓ Between-dataset Ensemble result saved to %s", basename(output_file)))

  return(output_file)
}
