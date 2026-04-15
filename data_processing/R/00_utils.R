# data_processing/R/00_utils.R --- Data processing utility functions
#
# NOTE: This file is sourced by _targets.R which already loads shared_helpers.R and cli_utils.R
# No need to source those files again here; they'll be in the global environment

# ============================================================================
# LOAD DATA PROCESSING PARAMETERS
# ============================================================================

load_dp_params <- function() {
  # Load parameters from data_processing/config/processing_parameters.yaml
  # Use direct path since load_pipeline_config() doesn't support data_processing module
  root <- proj_root()
  config_path <- file.path(root, "data_processing", "config", "processing_parameters.yaml")
  
  if (!file.exists(config_path)) {
    stop("Config file not found: ", config_path)
  }
  
  config <- yaml::read_yaml(config_path)

  # Make paths robust to current working directory by resolving
  # them relative to the project root.
  if (!is.null(config$in_dir)) {
    config$in_dir <- proj_path(config$in_dir)
  }
  if (!is.null(config$out_dir)) {
    config$out_dir <- proj_path(config$out_dir)
  }
  
  # Ensure output directory exists
  dir.create(config$out_dir, showWarnings = FALSE, recursive = TRUE)
  config
}

load_dp_datasets <- function() {
  # Load dataset registry from data_processing/config/core_datasets.yaml
  # Use direct path since load_pipeline_config() doesn't support data_processing module
  root <- proj_root()
  config_path <- file.path(root, "data_processing", "config", "core_datasets.yaml")
  
  if (!file.exists(config_path)) {
    stop("Config file not found: ", config_path)
  }
  
  `%||%` <- function(x, y) if (is.null(x)) y else x

  datasets_config <- yaml::read_yaml(config_path)$datasets
  if (is.null(datasets_config) || length(datasets_config) == 0) {
    stop("No datasets found in: ", config_path)
  }

  # Convert YAML list entries into a data.frame for targets::pattern = map().
  # Keep list-like fields as list-columns.
  df <- data.frame(
    id = vapply(datasets_config, function(x) x$id %||% NA_character_, character(1)),
    reference = vapply(datasets_config, function(x) x$reference %||% NA_character_, character(1)),
    ident_col = vapply(datasets_config, function(x) x$ident_col %||% NA_character_, character(1)),
    sample_col = vapply(datasets_config, function(x) x$sample_col %||% NA_character_, character(1)),
    batch_col = vapply(datasets_config, function(x) x$batch_col %||% "", character(1)),
    condition_col = vapply(datasets_config, function(x) x$condition_col %||% "", character(1)),
    stringsAsFactors = FALSE
  )

  df$raw_file <- I(lapply(datasets_config, function(x) x$raw_file %||% NULL))
  df$exclude_cell_types <- I(lapply(datasets_config, function(x) x$exclude_cell_types %||% character(0)))

  df
}

# ============================================================================
# DATA PROCESSING PATHS
# ============================================================================

dp_output_dir <- function(params = NULL) {
  if (is.null(params)) {
    params <- load_dp_params()
  }
  ensure_dir(proj_path(params$out_dir))
  proj_path(params$out_dir)
}

# ============================================================================
# UTILITY FUNCTIONS
# ============================================================================

#' Check if dataset output files exist
#' @param prefix Dataset prefix (e.g., "JoaI")
#' @param out_dir Output directory
#' @return TRUE if at least one output file exists, FALSE otherwise
output_exists <- function(prefix, out_dir) {
  files <- list.files(out_dir, pattern = paste0("^", prefix, "_.*\\.rds$"))
  length(files) > 0
}
