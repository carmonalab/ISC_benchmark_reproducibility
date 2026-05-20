suppressPackageStartupMessages({
  library(targets)
  library(tarchetypes)
  library(dplyr)
  library(tidyr)
})

source("../R/shared_helpers.R")
source(proj_path("R/cli_utils.R"))
source("R/resources_utils.R")

tar_option_set(
  packages = c(
    "yaml", "dplyr", "tidyr", "Matrix",
    "scTypeEval", "Seurat", "SeuratObject", "bench"
  ),
  storage = "main",
  retrieval = "main",
  error = "continue"
)

list(
  tar_target(
    resource_params,
    load_resource_config()
  ),

  tar_target(
    resource_requested_dataset_ids,
    get_requested_resource_dataset_ids()
  ),

  tar_target(
    resource_dataset_info,
    get_resource_dataset_idents(
      config = resource_params,
      selected_dataset_ids = resource_requested_dataset_ids
    )
  ),

  tar_target(
    resource_ident_grid,
    {
      ident_grid <- build_resource_ident_grid(resource_dataset_info)
      # Filter out already completed combinations
      filter_incomplete_ident_grid(ident_grid, resource_params)
    }
  ),

  tar_target(
    resource_metric_grid,
    build_resource_metric_grid(resource_params)
  ),

  tar_target(
    resource_prepared_input,
    prepare_resource_input(
      dataset_id = resource_ident_grid$dataset_id,
      dataset_file = resource_ident_grid$dataset_file,
      ident = resource_ident_grid$ident,
      params = resource_params
    ),
    format = "file",
    pattern = map(resource_ident_grid),
    iteration = "list"
  ),

  tar_target(
    resource_metric_result,
    benchmark_resource_pair(
      prepared_path = resource_prepared_input,
      dissimilarity_method = resource_metric_grid$dissimilarity_method,
      consistency_metric = resource_metric_grid$consistency_metric,
      params = resource_params
    ),
    format = "file",
    pattern = cross(resource_prepared_input, resource_metric_grid),
    iteration = "list"
  )
)