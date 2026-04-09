# ISC_benchmark/R/00_utils.R --- ISC-specific utility functions

source("../R/shared_helpers.R")
source("../R/cli_utils.R")

isc_params_path <- "config/isc_benchmark_parameters.yaml"

# ============================================================================
# LOAD ISC PARAMETERS
# ============================================================================

load_isc_params <- function(path = isc_params_path) {
  # Load all parameters from ISC_benchmark/config/isc_benchmark_parameters.yaml
  # (All global and ISC-specific settings are here)
  load_pipeline_config(path)
}

# ============================================================================
# ISC-SPECIFIC PATHS (from configuration)
# ============================================================================

isc_data_processed_dir <- function(params = NULL) {
  if (is.null(params)) {
    params <- load_isc_params()
  }
  proj_path(params$paths$data_processed)
}

isc_results_root <- function(params = NULL) {
  if (is.null(params)) {
    params <- load_isc_params()
  }
  root <- proj_path(params$paths$results_root)
  ensure_dir(root)
  root
}

get_isc_tasks <- function(params = NULL) {
  if (is.null(params)) {
    params <- load_isc_params()
  }
  # Return names of tasks to run (those with TRUE in params$isc_tasks$run)
  tasks <- params$run
  names(tasks)[which(unlist(tasks) == TRUE)]
}

# ============================================================================
# VALIDATION
# ============================================================================

validate_processed_isc_data <- function(dataset_id) {
  data_path <- file.path(isc_data_processed_dir(), paste0(dataset_id, ".rds"))
  meta_path <- paste0(data_path, ".yaml")
  
  if (!file.exists(data_path)) {
    warning("Processed ISC data missing: ", dataset_id)
    return(FALSE)
  }
  
  TRUE
}
