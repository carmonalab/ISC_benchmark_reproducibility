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
  
  datasets_config <- yaml::read_yaml(config_path)
  datasets_config$datasets
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
