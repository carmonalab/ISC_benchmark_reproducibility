#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript scripts/run_group_size.R <group_size>")
}

group_size <- as.integer(args[[1]])
if (is.na(group_size) || group_size < 1) {
  stop("group_size must be a positive integer")
}

source("R/00_utils.R")
source("R/saturation_helpers.R")

cfg <- load_saturation_config("config/saturation_parameters.yaml")
merged_input <- load_merged_isc_input(cfg)
saturation_input <- prepare_saturation_input(merged_input, cfg)

res <- run_saturation_group_size(group_size, saturation_input, cfg)
cat(sprintf("Completed group size %d. Cached rows: %d\n", group_size, nrow(res)))
