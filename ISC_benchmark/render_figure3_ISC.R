#!/usr/bin/env Rscript

stop_if_not_root <- function() {
  if (!file.exists("config/benchmark_parameters.yaml") || !dir.exists("ISC_benchmark")) {
    stop("Run this script from the repository root.")
  }
}

main <- function() {
  stop_if_not_root()

  if (!requireNamespace("rmarkdown", quietly = TRUE)) {
    stop("Package 'rmarkdown' is required to render the figure.")
  }

  out_dir <- file.path("results", "figures")
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  input <- file.path("ISC_benchmark", "figure3_ISC.Rmd")
  if (!file.exists(input)) stop("Missing input Rmd: ", input)

  rmarkdown::render(
    input = input,
    output_dir = out_dir,
    output_file = "figure3_ISC.html",
    knit_root_dir = normalizePath(getwd(), winslash = "/", mustWork = TRUE),
    quiet = FALSE
  )
}

if (!interactive()) {
  main()
}
