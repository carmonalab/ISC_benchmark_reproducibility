# ISC_benchmark/_targets.R
#
# HPC-optimized targets pipeline for ISC benchmark
#
# Features:
#   - Separate target for each dataset+ident+task combination
#   - Designed for HPC job submission: each target = one job
#   - Failure tracking and selective retries capability
#   - Incremental processing: only reruns missing/failed outputs
#   - Error isolation: one task failure doesn't halt pipeline
#
# Task progression:
#   1-2: Sensitivity to signal degradation and over-partitioning
#   3-6: Robustness to annotation granularity, complexity, dataset size
#   7-8: Robustness to batch effects and biological perturbations
#
# Architecture:
#   - Static targets: config, dataset list, active tasks
#   - Dynamic targets: expand_grid(datasets × tasks) with nested ident tar_map
#   - Result: ~100+ independent targets (fine-grained parallelization)
#
# Usage:
#   cd ISC_benchmark
#   targets::tar_make()                # Run all pending targets
#   targets::tar_status()              # Check which tasks succeeded/failed
#   targets::tar_make()                # Rerun only failed tasks (automatic)
#
# HPC Submission:
#   bash scripts/submit_hpc.sh         # Submits one SLURM job per dataset
#
# NOTE: Run from ISC_benchmark directory

library(targets)
library(yaml)
library(dplyr)
library(tidyr)

# Source utilities and task functions
# (project_paths.R deleted — all functions moved to shared_helpers.R)
source("../R/shared_helpers.R")
source(proj_path("R/cli_utils.R"))
source("R/isc_benchmark_helpers.R")
source("R/isc_tasks.R")

# ============================================================================
# TARGETS PIPELINE CONFIGURATION
# ============================================================================

tar_option_set(
  packages = c(
    "yaml", "dplyr", "tidyr", "Seurat", "SeuratObject", 
    "Matrix", "scTypeEval", "BiocParallel"
  ),
  storage = "worker",
  retrieval = "worker"
)

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

#' Load ISC benchmark configuration from YAML
#'
#' @return List with configuration, paths, and options
get_isc_config <- function() {
  config <- read_yaml("config/isc_benchmark_parameters.yaml")
  
  # Set absolute paths (relative to ISC_benchmark directory)
  config$processed_data_dir <- normalizePath("../data/processed")
  config$output$dir <- normalizePath("../results/isc_benchmark")
  config$dataset_idents_file <- normalizePath("config/dataset_idents.yaml")
  
  # Create output directories
  dir.create(config$output$dir, showWarnings = FALSE, recursive = TRUE)
  
  config
}

#' Get datasets and their cell type annotation columns
#'
#' Reads dataset_idents.yaml and matches with processed data stems.
#' Returns a table suitable for tar_map with dataset_id, dataset_file, ident_cols.
#'
#' @param config Configuration list
#' @param selected_dataset_ids Optional character vector of processed dataset stems to keep
#' @param selected_dataset_families Optional character vector of dataset families to keep
#' @return Data frame with dataset_id, dataset_file, dataset_family, dataset_stems, ident_cols
normalize_dataset_family <- function(dataset_id) {
  family <- sub("_.*$", "", dataset_id)
  if (family == "StephensonE") {
    return("Stephenson")
  }
  family
}

get_dataset_idents <- function(config, selected_dataset_ids = NULL, selected_dataset_families = NULL) {
  ident_mapping <- read_yaml(config$dataset_idents_file)$idents

  proc_dir <- config$processed_data_dir
  if (!dir.exists(proc_dir)) {
    stop("Processed data directory not found: ", proc_dir)
  }

  files <- list.files(proc_dir, pattern = "\\.rds$", full.names = TRUE)
  if (length(files) == 0) {
    stop("No processed datasets found in: ", proc_dir)
  }

  basename_only <- basename(files)
  dataset_ids <- unique(sub("\\.rds$", "", basename_only))

  dataset_info <- lapply(dataset_ids, function(dataset_id) {
    dataset_file <- file.path(proc_dir, paste0(dataset_id, ".rds"))

    if (!file.exists(dataset_file)) {
      warning("No file found for dataset: ", dataset_id)
      return(NULL)
    }

    dataset_family <- normalize_dataset_family(dataset_id)
    ident_cols <- ident_mapping[[dataset_family]]

    if (is.null(ident_cols)) {
      warning("No ident columns defined in config for: ", dataset_family)
      return(NULL)
    }

    data.frame(
      dataset_id = dataset_id,
      dataset_file = dataset_file,
      dataset_family = dataset_family,
      dataset_stems = dataset_id,
      ident_cols = paste(ident_cols, collapse = ","),
      stringsAsFactors = FALSE
    )
  })

  dataset_info <- Filter(Negate(is.null), dataset_info)
  if (length(dataset_info) == 0) {
    stop("No valid datasets found after checking configurations")
  }

  dataset_info <- do.call(rbind, dataset_info)

  if (!is.null(selected_dataset_ids)) {
    selected_dataset_ids <- unique(trimws(selected_dataset_ids))
    selected_dataset_ids <- selected_dataset_ids[nzchar(selected_dataset_ids)]

    if (length(selected_dataset_ids) > 0) {
      missing_ids <- setdiff(selected_dataset_ids, dataset_info$dataset_id)
      if (length(missing_ids) > 0) {
        stop(
          "Requested dataset(s) not available in processed inputs: ",
          paste(missing_ids, collapse = ", ")
        )
      }

      dataset_info <- dataset_info[match(selected_dataset_ids, dataset_info$dataset_id), , drop = FALSE]
    }
  }

  if (!is.null(selected_dataset_families)) {
    selected_dataset_families <- unique(trimws(selected_dataset_families))
    selected_dataset_families <- selected_dataset_families[nzchar(selected_dataset_families)]

    if (length(selected_dataset_families) > 0) {
      missing_families <- setdiff(selected_dataset_families, unique(dataset_info$dataset_family))
      if (length(missing_families) > 0) {
        stop(
          "Requested dataset family(ies) not available in processed inputs: ",
          paste(missing_families, collapse = ", ")
        )
      }

      family_rows <- lapply(selected_dataset_families, function(family_id) {
        family_matches <- dataset_info[dataset_info$dataset_family == family_id, , drop = FALSE]
        if (nrow(family_matches) == 0) {
          return(NULL)
        }

        data.frame(
          dataset_id = family_id,
          dataset_file = family_matches$dataset_file[1],
          dataset_family = family_id,
          dataset_stems = paste(family_matches$dataset_id, collapse = ","),
          ident_cols = family_matches$ident_cols[1],
          stringsAsFactors = FALSE
        )
      })

      family_rows <- Filter(Negate(is.null), family_rows)
      if (length(family_rows) > 0) {
        dataset_info <- do.call(rbind, family_rows)
      } else {
        dataset_info <- dataset_info[0, , drop = FALSE]
      }
    }
  }

  dataset_info
}

#' Read requested dataset IDs from the environment
#'
#' Supports comma-separated selections via ISC_DATASET_IDS or single-dataset
#' test runs via ISC_DATASET_ID / ISC_TEST_DATASET.
#'
#' @return Character vector of requested dataset IDs, or NULL for all datasets
get_requested_dataset_ids <- function() {
  requested <- c(
    Sys.getenv("ISC_DATASET_IDS", unset = ""),
    Sys.getenv("ISC_DATASET_ID", unset = ""),
    Sys.getenv("ISC_TEST_DATASET", unset = "")
  )

  requested <- requested[nzchar(requested)]
  if (length(requested) == 0) {
    return(NULL)
  }

  unique(trimws(unlist(strsplit(paste(requested, collapse = ","), ",", fixed = TRUE))))
}

#' Read requested dataset families from the environment
#'
#' Supports comma-separated selections via ISC_DATASET_FAMILIES or ISC_DATASET_FAMILY.
#'
#' @return Character vector of requested dataset families, or NULL for all families
get_requested_dataset_families <- function() {
  requested <- c(
    Sys.getenv("ISC_DATASET_FAMILIES", unset = ""),
    Sys.getenv("ISC_DATASET_FAMILY", unset = "")
  )

  requested <- requested[nzchar(requested)]
  if (length(requested) == 0) {
    return(NULL)
  }

  unique(trimws(unlist(strsplit(paste(requested, collapse = ","), ",", fixed = TRUE))))
}

#' Read requested task names from the environment
#'
#' Supports comma-separated task lists via ISC_TASKS.
#'
#' @return Character vector of requested tasks, or NULL for all enabled tasks
get_requested_task_names <- function() {
  requested <- Sys.getenv("ISC_TASKS", unset = "")
  requested <- trimws(requested)
  if (!nzchar(requested)) {
    return(NULL)
  }

  unique(trimws(unlist(strsplit(requested, ",", fixed = TRUE))))
}
#' Get active task names from config
#'
#' @param config Configuration list
#' @return Character vector of task names
get_active_tasks <- function(config) {
  tasks <- names(config$run)[config$run == TRUE]
  requested_tasks <- get_requested_task_names()

  if (!is.null(requested_tasks)) {
    tasks <- intersect(tasks, requested_tasks)
  }

  if (length(tasks) == 0) {
    stop("No enabled tasks selected for execution")
  }

  tasks
}

#' Parse comma-separated ident columns string
#'
#' @param ident_cols_str Comma-separated string (e.g., "cell.type,annotation")
#' @return Character vector
parse_ident_cols <- function(ident_cols_str) {
  trimws(strsplit(ident_cols_str, ",")[[1]])
}

# ============================================================================
# TARGETS PIPELINE (HPC-OPTIMIZED)
# ============================================================================
# Each dataset+ident+task = separate target = separate HPC job
# Enables fine-grained parallelization and failure tracking

list(
  # ========== STATIC CONFIGURATION TARGETS ==========
  
  # Load configuration (tracks file changes)
  tar_target(
    config,
    get_isc_config(),
    cue = tar_cue(file = "config/isc_benchmark_parameters.yaml")
  ),
  
  # Get dataset-ident mappings (tracks dataset_idents.yaml changes)
  tar_target(
    datasets_idents,
    get_dataset_idents(
      config,
      selected_dataset_ids = get_requested_dataset_ids(),
      selected_dataset_families = get_requested_dataset_families()
    ),
    cue = tar_cue(file = "config/dataset_idents.yaml")
  ),
  
  # Get active task names
  tar_target(
    active_tasks,
    get_active_tasks(config)
  ),
  
  # ========== DYNAMIC TARGETS: Per dataset × task × ident ==========
  # Each target is one independent task execution
  # expand_grid creates all combinations, nested tar_map handles ident arrays
  
  tar_map(
    values = expand_grid(datasets_idents, task_name = active_tasks),
    
    # Parse ident columns string into vector
    tar_target(
      name = idents_vector,
      command = parse_ident_cols(ident_cols)
    ),
    
    # Create per-ident targets (one per cell type annotation)
    tar_map(
      values = list(ident_col = idents_vector),
      
      # Individual task execution target
      # Each represents one complete unit of work suitable for HPC job submission
      tar_target(
        name = task_result,
        command = {
          # Logging
          timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
          job_id <- sprintf("%s:%s:%s", dataset_id, ident_col, task_name)
          cat(sprintf("[%s] Starting job %s\n", timestamp, job_id))
          
          # Create output directory
          output_dir <- file.path(config$output$dir, dataset_id)
          dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
          
          # Execute single task
          result <- tryCatch({
            run_isc_benchmark_on_dataset(
              dataset_id = dataset_id,
              ident_col = ident_col,
              task_name = task_name,
              dataset_path = dataset_file,
              dataset_stems = dataset_stems,
              config = config,
              output_dir = output_dir
            )
          }, error = function(e) {
            # Return error data frame on failure
            data.frame(
              dataset_id = dataset_id,
              ident = ident_col,
              task = task_name,
              status = "failed",
              error = as.character(e$message),
              n_results = 0,
              stringsAsFactors = FALSE
            )
          })
          
          # Log completion
          cat(sprintf("[%s] Completed job %s - status: %s\n",
                      timestamp, job_id, result$status[1]))
          
          result
        },
        error = "continue"  # Continue pipeline even if one task fails
      )
    )
  ),
  
  # ========== AGGREGATION TARGETS ==========
  
  # Aggregate all task results into single data frame
  tar_target(
    all_results,
    {
      # Combine all task_result targets (targets framework handles pattern matching)
      all_task_results <- tar_read_raw_list(
        grep("^task_result_", tar_progress()$name, value = TRUE)
      )
      
      # Bind into single data frame
      combined <- do.call(rbind, all_task_results)
      rownames(combined) <- NULL
      
      message(sprintf("\nAggregated %d task results", nrow(combined)))
      combined
    },
    pattern = cross(grep("^task_result_", names(.targets), value = TRUE))
  ),
  
  # Save aggregated results
  tar_target(
    save_aggregated,
    {
      output_file <- file.path(config$output$dir, "all_results.rds")
      saveRDS(all_results, output_file)
      
      message("\nSaved aggregated results to: ", output_file)
      invisible(output_file)
    }
  ),
  
  # ========== FAILURE REPORTING ==========
  
  # Identify and report failed tasks
  tar_target(
    failure_report,
    {
      failures <- all_results[all_results$status == "failed", ]
      
      cat("\n", strrep("=", 70), "\n")
      cat("ISC BENCHMARK STATUS REPORT\n")
      cat(strrep("=", 70), "\n\n")
      
      if (nrow(failures) > 0) {
        cat("⚠ FAILED TASKS (", nrow(failures), "):\n\n")
        print(failures[, c("dataset_id", "ident", "task", "error")])
        
        # Save failure list for potential retry
        failure_file <- file.path(config$output$dir, "failures.csv")
        write.csv(failures, failure_file, row.names = FALSE)
        cat("\nFailure report saved to: ", failure_file, "\n")
        cat("To retry failed tasks, run: targets::tar_make() again\n")
      } else {
        cat("✓ ALL TASKS COMPLETED SUCCESSFULLY!\n")
      }
      
      # Summary by task
      cat("\n--- By Task ---\n")
      task_summary <- all_results %>%
        group_by(task, status) %>%
        summarise(n = n(), .groups = "drop") %>%
        pivot_wider(names_from = status, values_from = n, values_fill = 0)
      print(as.data.frame(task_summary))
      
      # Summary by dataset
      cat("\n--- By Dataset ---\n")
      dataset_summary <- all_results %>%
        group_by(dataset_id, status) %>%
        summarise(n = n(), .groups = "drop") %>%
        pivot_wider(names_from = status, values_from = n, values_fill = 0)
      print(as.data.frame(dataset_summary))
      
      cat("\nResults directory: ", config$output$dir, "\n")
      cat(strrep("=", 70), "\n\n")
      
      invisible(failures)
    }
  )
)
