#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(yaml)
})

source("utils/cli_utils.R")

`%||%` <- function(a, b) if (!is.null(a)) a else b

args <- parse_args(commandArgs(trailingOnly = TRUE))

processed_dir <- args[["processed_dir"]] %||% file.path("data", "processed", "isc")
params_path <- args[["params"]] %||% file.path("ISC_benchmark", "BatchEffect_singleRun", "params.yaml")
ncores <- as.numeric(args[["ncores"]] %||% 4)
dry_run <- tolower(args[["dry_run"]] %||% "no") %in% c("1", "true", "t", "yes", "y")

out_task7 <- args[["out_task7"]] %||% file.path("results", "isc", "tasks", "task7_batch_comparison", "output")
out_task8 <- args[["out_task8"]] %||% file.path("results", "isc", "tasks", "task8_perturbation_comparison", "output")

if (!dir.exists(processed_dir)) {
  stop("Processed dir not found: ", processed_dir)
}
if (!file.exists(params_path)) {
  stop("Params file not found: ", params_path)
}

rds_files <- list.files(processed_dir, pattern = "\\.rds$", full.names = TRUE, recursive = TRUE)
if (length(rds_files) == 0) {
  stop("No processed .rds files found under: ", processed_dir)
}

run_one <- function(dataset_path, split_col, out_dir, idents, label) {
  cmd_args <- c(
    file.path("ISC_benchmark", "BatchEffect_singleRun", "BatchEffect_singleRun.R"),
    "--dataset", dataset_path,
    "--params", params_path,
    "--idents", paste(idents, collapse = ","),
    "--ncores", as.character(ncores),
    "--split_col", split_col,
    "--out", out_dir
  )

  message_time(paste0("[", label, "] ", basename(dataset_path), " (split_col=", split_col, ")"))
  if (dry_run) {
    message("DRY RUN: Rscript ", paste(cmd_args, collapse = " "))
    return(invisible(NULL))
  }

  status <- system2("Rscript", args = cmd_args)
  if (!identical(status, 0L)) {
    stop("BatchEffect_singleRun failed (exit code ", status, ") for ", dataset_path)
  }
}

for (dataset_path in rds_files) {
  meta_path <- paste0(dataset_path, ".yaml")
  if (!file.exists(meta_path)) {
    message("Skipping (missing sidecar YAML): ", dataset_path)
    next
  }

  meta <- yaml::read_yaml(meta_path)
  idents <- meta$idents %||% character()
  idents <- idents[!is.na(idents)]
  idents <- idents[idents != ""]

  if (length(idents) == 0) {
    message("Skipping (no idents): ", dataset_path)
    next
  }

  if (isTRUE(meta$batch_comparison)) {
    run_one(dataset_path, split_col = "split", out_dir = out_task7, idents = idents, label = "Task 7")
  }

  if (isTRUE(meta$perturbation_comparison)) {
    run_one(dataset_path, split_col = "condition", out_dir = out_task8, idents = idents, label = "Task 8")
  }
}

message_time("Tasks 7–8 completed")
