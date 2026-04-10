#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(yaml)
})

load_saturation_config <- function(config_path = "config/saturation_parameters.yaml") {
  cfg <- yaml::read_yaml(config_path)

  cfg$paths$merged_results_file <- normalizePath(
    cfg$paths$merged_results_file,
    winslash = "/",
    mustWork = FALSE
  )
  cfg$paths$fallback_metrics_dir <- normalizePath(
    cfg$paths$fallback_metrics_dir,
    winslash = "/",
    mustWork = FALSE
  )
  cfg$paths$output_dir <- normalizePath(
    cfg$paths$output_dir,
    winslash = "/",
    mustWork = FALSE
  )

  dir.create(cfg$paths$output_dir, recursive = TRUE, showWarnings = FALSE)

  cfg
}
