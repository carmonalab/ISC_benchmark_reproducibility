#!/usr/bin/env Rscript

# Source shared utilities first to get find_proj_root() function
source("../R/shared_helpers.R")
project_root <- find_proj_root()

# Set up renv library paths
r_mm <- paste0(R.version$major, ".", sub("\\..*$", "", R.version$minor))
renv_lib <- file.path(project_root, "renv", "library", paste0("R-", r_mm), R.version$platform)
if (dir.exists(renv_lib)) {
  .libPaths(unique(c(renv_lib, .libPaths())))
}

suppressPackageStartupMessages({
  library(data.table)
})

source(proj_path("ISC_benchmark_saturation/R/00_utils.R"))
source(proj_path("ISC_benchmark_saturation/R/saturation_helpers.R"))

cfg <- load_saturation_config("config/saturation_parameters.yaml")
merged_input <- load_merged_isc_input(cfg)
saturation_input <- prepare_saturation_input(merged_input, cfg)

all_group <- lapply(cfg$analysis$group_sizes, function(k) {
  load_group_cache(k, cfg)
})

all_group <- all_group[vapply(all_group, function(x) is.data.frame(x) && nrow(x) > 0, logical(1))]
if (length(all_group) == 0) {
  stop("No cached group-size results found. Run group jobs first.")
}

combined <- data.table::rbindlist(all_group, fill = TRUE)
if ("ComparisonID" %in% names(combined)) {
  combined <- combined[!duplicated(combined$ComparisonID), ]
}

out <- save_saturation_outputs(saturation_input, combined, cfg)
cat(sprintf("Saved saturation outputs:\n- %s\n- %s\n", out$input_file, out$results_file))
