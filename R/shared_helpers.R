# shared_helpers.R --- Utilities shared by both ISC_benchmark and label_transfer_task

library(tidyverse)
library(yaml)

# ============================================================================
# PATH HELPERS (root-agnostic: work from any subdirectory)
# ============================================================================

find_proj_root <- function(start_dir = getwd()) {
  # Walk up directory tree until we find project root marker
  # (look for .Rproj file or renv.lock at root)
  # Each pipeline subdirectory has its own config/, but projects root is where
  # data/ and renv.lock live.
  
  current <- normalizePath(start_dir, mustWork = TRUE)
  max_depth <- 10
  depth <- 0
  
  while (depth < max_depth) {
    # Look for .Rproj file or renv.lock
    rproj_files <- list.files(current, pattern = "\\.Rproj$")
    if (length(rproj_files) > 0 || file.exists(file.path(current, "renv.lock"))) {
      return(current)
    }
    parent <- dirname(current)
    if (parent == current) break  # Hit filesystem root
    current <- parent
    depth <- depth + 1
  }
  
  stop("Could not find project root (expected .Rproj file or renv.lock). ",
       "Are you in ISC_benchmark_reproducibility/ISC_benchmark/ or ",
       "ISC_benchmark_reproducibility/label_transfer_task/?")
}

proj_root <- function() {
  find_proj_root()
}

proj_path <- function(...) {
  file.path(proj_root(), ...)
}

ensure_dir <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
  invisible(path)
}

# ============================================================================
# GET CURRENT PIPELINE (returns "ISC_benchmark" or "label_transfer_task")
# ============================================================================

get_current_pipeline <- function() {
  # Determine which pipeline is running based on working directory
  cwd <- getwd()
  
  if (grepl("ISC_benchmark/?$", cwd)) {
    return("ISC_benchmark")
  } else if (grepl("label_transfer_task/?$", cwd)) {
    return("label_transfer_task")
  } else {
    stop("Unable to determine current pipeline. ",
         "Expected working directory to be ISC_benchmark/ or label_transfer_task/")
  }
}

# ============================================================================
# CONFIGURATION LOADING
# ============================================================================

load_config <- function(path) {
  if (!file.exists(path)) {
    stop("Config file not found: ", path)
  }
  yaml::read_yaml(path)
}

load_pipeline_config <- function(filename) {
  # Load from current pipeline's config/ subdirectory
  # All required parameters should be in pipeline-specific config
  pipeline <- get_current_pipeline()
  load_config(proj_path(pipeline, "config", filename))
}

load_data_processing_config <- function(filename) {
  # Load from data_processing/config/ subdirectory
  load_config(proj_path("data_processing", "config", filename))
}

load_dataset_registry <- function() {
  # Load datasets from data_processing/config/datasets.yaml
  load_data_processing_config("datasets.yaml")$datasets
}

# ============================================================================
# REPRODUCIBILITY HELPERS
# ============================================================================

compute_file_hash <- function(path) {
  if (!file.exists(path)) return(NA_character_)
  tools::md5sum(path) %>% as.character() %>% unname()
}
