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
    ncores = 1) {

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
    warning("Classifier failed: ", e$message)
    NA
  })

  # Ensure predictions can be aligned downstream
  if (!all(is.na(predictions)) && length(predictions) == ncol(query$counts)) {
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
