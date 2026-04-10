# label_transfer_task/R/00_utils.R --- Label-transfer-specific utility functions

source("../R/shared_helpers.R")

# ============================================================================
# LABEL-TRANSFER-SPECIFIC PATHS
# ============================================================================

lt_data_processed_dir <- function() {
  proj_path("data/processed/label_transfer")
}

lt_results_root <- function() {
  ensure_dir(proj_path("label_transfer_task/results"))
  proj_path("label_transfer_task/results")
}

lt_raw_results_dir <- function() {
  ensure_dir(file.path(lt_results_root(), "raw_results"))
}

lt_aggregated_dir <- function() {
  ensure_dir(file.path(lt_results_root(), "aggregated"))
}

lt_figures_dir <- function() {
  ensure_dir(file.path(lt_results_root(), "figures"))
}

# ============================================================================
# LOAD LABEL-TRANSFER PARAMETERS
# ============================================================================

load_lt_params <- function() {
  # Load all label transfer parameters from pipeline-specific config
  # (All required settings are self-contained in label_transfer_parameters.yaml)
  load_pipeline_config("label_transfer_parameters.yaml")
}

get_lt_classifiers <- function(params) {
  # Return list of classifiers to evaluate
  params$classifiers$methods %||% c("SingleR", "RandomForest", "SVM")
}

# ============================================================================
# DATASET IDENTIFICATION COLUMN MAPPING
# ============================================================================

#' Get cell-type identification column name by dataset prefix
#'
#' Maps dataset identifiers to their corresponding cell-type annotation columns.
#' This ensures consistent identification across different datasets.
#'
#' @param prefix Character: Dataset prefix (e.g., "JoaI", "Mitchel", "BCC")
#'
#' @return Character: Column name to use for cell-type identification
#'
#' @details
#' Maps dataset prefixes to their metadata columns:
#' - JoaI (Joanito et al., 2022 CRC): "cell.type"
#' - Mitchel (Mitchel et al., 2023): "OriginalAnnotationLevel1"
#' - BCC (Yost et al., BCC atlas): "annotation"
#' - LungAtlas (Sikkema et al., HCA Lung): "cell_type"
#' - ICBAtlas (Gondal et al., ICB Atlas): "cell_type"
#'
get_idents_by_prefix <- function(prefix) {
  switch(prefix,
    "JoaI"        = "cell.type",
    "Mitchel"     = "OriginalAnnotationLevel1",
    "Joanito"     = "cell.type",
    "Lee"         = "cell.type",
    "Stephenson"  = "OriginalAnnotationLevel2",
    "BCC"         = "annotation",
    "LungAtlas"   = "cell_type",
    "ICBAtlas"    = "cell_type",
    "Yerly"       = "annotation",
    "Ganier"      = "annotation",
    "celltype"    # Default fallback
  )
}

# ============================================================================
# VALIDATION & DATA PREPARATION
# ============================================================================

validate_label_transfer_participation <- function(dataset_id) {
  # Check if dataset participates in label_transfer task
  meta_path <- file.path(
    proj_path("data/processed/isc"),
    paste0(dataset_id, ".rds.yaml")
  )
  
  if (!file.exists(meta_path)) {
    return(TRUE)  # Default: include if metadata missing
  }
  
  meta <- load_config(meta_path)
  isTRUE(meta$label_transfer)
}

validate_label_transfer_data <- function(dataset_id) {
  # Check if query + reference files exist
  data_dir <- lt_data_processed_dir()
  query_path <- file.path(data_dir, dataset_id, "query.rds")
  ref_path <- file.path(data_dir, dataset_id, "reference.rds")
  
  all(file.exists(query_path, ref_path))
}

# ============================================================================
# HELPER: LOAD CLASSIFIER DEFINITIONS
# ============================================================================

load_classifier_code <- function() {
  # Source the classifiers.R file that defines classifier functions
  classifier_file <- proj_path("label_transfer_task/classifiers/classifiers.R")
  
  if (!file.exists(classifier_file)) {
    stop("Classifier definitions not found: ", classifier_file)
  }
  
  source(classifier_file, local = TRUE)
  invisible(NULL)
}
