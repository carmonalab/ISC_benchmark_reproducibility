#!/usr/bin/env Rscript

library(targets)
library(tibble)
library(data.table)

# Source shared utilities (defines proj_path)
source("../R/shared_helpers.R")
source(proj_path("ISC_benchmark_saturation/R/00_utils.R"))
source(proj_path("ISC_benchmark_saturation/R/saturation_helpers.R"))

tar_option_set(
  packages = c("yaml", "dplyr", "tibble", "tidyr", "purrr", "stringr", "data.table", "BiocParallel")
)

list(
  tar_target(
    cfg,
    load_saturation_config("config/saturation_parameters.yaml"),
    cue = tar_cue(file = "config/saturation_parameters.yaml")
  ),

  tar_target(
    merged_input,
    load_merged_isc_input(cfg)
  ),

  tar_target(
    saturation_input,
    prepare_saturation_input(merged_input, cfg)
  ),

  tar_target(
    group_size,
    cfg$analysis$group_sizes
  ),

  tar_target(
    saturation_by_group,
    run_saturation_group_size(group_size, saturation_input, cfg),
    pattern = map(group_size)
  ),

  tar_target(
    saturation_results,
    data.table::rbindlist(saturation_by_group, fill = TRUE)
  ),

  tar_target(
    save_outputs,
    save_saturation_outputs(saturation_input, saturation_results, cfg)
  )
)
