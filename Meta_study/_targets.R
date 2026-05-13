# Meta_study/_targets.R --- Targets pipeline for Meta-study consistency objects

suppressPackageStartupMessages({
  library(targets)
  library(tarchetypes)
  library(yaml)
  library(dplyr)
  library(tibble)
})

# Shared root/path helpers
source("../R/shared_helpers.R")

# Pipeline-specific helper functions
source("R/00_utils.R")

tar_option_set(
  packages = c("yaml", "dplyr", "tibble", "Matrix"),
  storage = "main",
  retrieval = "main"
)

list(
  tar_target(
    meta_params_file,
    "config/meta_study_parameters.yaml",
    format = "file"
  ),

  tar_target(
    meta_datasets_file,
    "config/datasets_metadata.yaml",
    format = "file"
  ),

  tar_target(
    meta_params,
    {
      meta_params_file
      load_meta_params()
    }
  ),

  tar_target(
    meta_datasets,
    {
      meta_datasets_file
      load_meta_datasets()
    }
  ),

  tar_target(
    meta_api_ready,
    {
      ensure_scTypeEval_api(meta_params)
      TRUE
    }
  ),

  tar_target(
    meta_black_list,
    {
      meta_api_ready
      load_black_list_vector(meta_params)
    }
  ),

  tar_target(
    meta_dissimilarity_methods,
    get_meta_dissimilarity_methods(meta_params)
  ),

  tar_target(
    meta_jobs,
    build_meta_jobs(meta_datasets)
  ),

  tar_target(
    meta_scTypeEval_file,
    {
      meta_api_ready
      run_meta_scTypeeval_job(
        ds_row = meta_jobs,
        ident_col = meta_jobs$ident_col,
        params = meta_params,
        black_list = meta_black_list,
        dissimilarity_methods = meta_dissimilarity_methods
      )
    },
    pattern = map(meta_jobs),
    iteration = "list",
    format = "file",
    error = "continue"
  ),

  tar_target(
    meta_summary_file,
    write_meta_summary(meta_scTypeEval_file, meta_params),
    format = "file"
  )
)
