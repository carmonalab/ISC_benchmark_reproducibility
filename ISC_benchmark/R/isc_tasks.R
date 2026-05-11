#!/usr/bin/env Rscript
# ISC_benchmark/R/01_isc_tasks.R
#
# ISC benchmark task definitions (reorganized by task progression)
#
# Tasks are ordered to progress from:
# (1-2) Sensitivity tests (robustness to signal degradation)
# (3-6) Robustness tests (annotation/complexity/dataset size)
# (7-8) Robustness tests (batch/biological effects)
#
# Each task calls scTypeEval wr_* functions and extracts metrics

source("R/isc_benchmark_helpers.R")

# ============================================================================
# TASK EXECUTION ORCHESTRATOR
# ============================================================================

#' Run ISC benchmark on single dataset and task
#'
#' Executes one task on one dataset with specified replicates and parameters
#'
#' @param dataset_id Dataset identifier (e.g., "JoaI")
#' @param ident_col Cell type annotation column name
#' @param task_name Task name from catalog (e.g., "missclassify")
#' @param dataset_path Path to processed dataset (Seurat object)
#' @param config Configuration list from YAML
#' @param output_dir Output directory for results
#'
#' @return Data frame with task results and metrics
run_isc_benchmark_on_dataset <- function(dataset_id,
                                         ident_col,
                                         task_name,
                                         dataset_path,
                                         dataset_stems = dataset_id,
                                         config,
                                         output_dir) {
  
  set.seed(config$seed)
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  # Use default blacklist (TCR, Ig, Y-genes) if none specified in config
  if (is.null(config$common$black_list)) {
    config$common$black_list <- get_default_blacklist()
  }

  
  message_step("Running ISC_BENCHMARK",
               sprintf("dataset=%s, task=%s, ident=%s", dataset_id, task_name, ident_col))
  
  # ========== STEP 1: Load dataset ==========
  obj <- tryCatch({
    if (task_name %in% c("batch_effects", "biological_perturbations")) {
      loaded_obj <- load_and_merge_processed_datasets(dataset_stems, config)
      message_step("LOAD", sprintf("Merged dataset family loaded: %d cells, %d genes from %d stem(s)",
                                    ncol(loaded_obj), nrow(loaded_obj),
                                    length(unique(trimws(unlist(strsplit(paste(dataset_stems, collapse = ","), ",", fixed = TRUE)))))))
      loaded_obj
    } else {
      if (!file.exists(dataset_path)) {
        stop("Dataset not found: ", dataset_path)
      }
      loaded_obj <- readRDS(dataset_path)
      message_step("LOAD", sprintf("Dataset loaded: %d cells, %d genes", ncol(loaded_obj), nrow(loaded_obj)))
      loaded_obj
    }
  }, error = function(e) {
    message_step("ERROR", sprintf("Dataset loading failed: %s", e$message))
    NULL
  })

  if (is.null(obj)) {
    return(data.frame(
      dataset_id = dataset_id,
      ident = ident_col,
      task = task_name,
      status = "failed",
      error = "load_failed",
      n_results = 0
    ))
  }
  
  # ========== STEP 2: Prepare for scTypeEval ==========
  obj_prepared <- tryCatch({
    prepare_scTypeEval_object(obj, ident_col, config)
  }, error = function(e) {
    message_step("ERROR", sprintf("Preparation failed: %s", e$message))
    return(NULL)
  })
  
  if (is.null(obj_prepared)) {
    return(data.frame(
      dataset_id = dataset_id,
      ident = ident_col,
      task = task_name,
      status = "failed",
      error = "preparation_failed",
      n_results = 0
    ))
  }
  
  # ========== STEP 3: Execute task ==========
  # ========== STEP 2b: Load/compute unified baseline dataframe for tasks 1-6 ==========
  # One full-dataset ISC computation (no perturbation) is shared across ALL tasks 1-6.
  # Cache and reuse only consistency dataframe (no full scTypeEval object persistence).
  TASKS_WITH_BASELINE <- c("missclassify", "SplitCelltype", "Nsamples", "NCell",
                            "Nct", "cellular_complexity")
  baseline_df  <- NULL

  if (task_name %in% TASKS_WITH_BASELINE) {
    baseline_cache_path <- file.path(output_dir, paste0("baseline_isc_", ident_col, ".rds"))

    baseline_df <- tryCatch(
      get_or_compute_baseline(obj_prepared, config, baseline_cache_path),
      error = function(e) {
        message_step("BASELINE", sprintf("Baseline computation failed (%s); task will recompute baseline internally", e$message))
        NULL
      }
    )
  }

  # ========== STEP 3: Execute task ==========
  task_output_dir <- file.path(output_dir, sprintf("%s_%s", task_name, ident_col))
  
  wr_result <- NULL
  task_metrics <- NULL
  
  tryCatch({
    message_step("TASK", sprintf("Running %s...", task_name))
    
    # Extract task-specific configuration from main config
    task_config_key <- paste0("task_", task_name)
    task_config <- config[[task_config_key]][["isc_params"]]
    metric_config <- config[[task_config_key]][["type"]]
    
    # Tasks 7/8 use the full sub-config directly and have no isc_params block;
    # only require isc_params for the tasks that actually use task_config.
    tasks_requiring_isc_params <- c("missclassify", "SplitCelltype", "Nct",
                                    "cellular_complexity", "Nsamples", "NCell")
    if (is.null(task_config) && task_name %in% tasks_requiring_isc_params) {
      stop("Task configuration not found: ", task_config_key)
    }
    
    # Dispatch to appropriate task function
    switch(task_name,
      "missclassify" = {
        wr_result <<- run_task_missclassify(obj_prepared, config, task_config, task_output_dir,
                                             baseline_df = baseline_df)
        task_metrics <<- extract_task_metrics(wr_result, task_name, metric_config)
      },
      "SplitCelltype" = {
        wr_result <<- run_task_SplitCelltype(obj_prepared, config, task_config, task_output_dir,
                                              baseline_df = baseline_df)
        task_metrics <<- extract_task_metrics(wr_result, task_name, metric_config)
      },
      "Nct" = {
        wr_result <<- run_task_Nct(obj_prepared, config, task_config, task_output_dir,
                                    baseline_df = baseline_df)
        task_metrics <<- extract_task_metrics(wr_result, task_name, metric_config)
      },
      "cellular_complexity" = {
        wr_result <<- run_task_cellular_complexity(obj_prepared, config, task_config, task_output_dir,
                                                   baseline_df = baseline_df)
        # Extract metrics from each complexity level
        task_metrics <<- extract_task_metrics(wr_result, task_name, metric_config)
      },
      "Nsamples" = {
        wr_result <<- run_task_Nsamples(obj_prepared, config, task_config, task_output_dir,
                                         baseline_df = baseline_df)
        task_metrics <<- extract_task_metrics(wr_result, task_name, metric_config)
      },
      "NCell" = {
        wr_result <<- run_task_NCell(obj_prepared, config, task_config, task_output_dir,
                                      baseline_df = baseline_df)
        task_metrics <<- extract_task_metrics(wr_result, task_name, metric_config)
      },
      "batch_effects" = {
        specs_file <- file.path(proj_root(), "data_processing", "config", "specs_datasets.csv")
        wr_result <<- run_task_batch_effects(obj_prepared, config, config[["task_batch_effects"]], task_output_dir,
                                             specs_path    = specs_file,
                                             results_root  = config$output$dir,
                                             dataset_stems = dataset_stems)
        if (!is.null(wr_result)) {
          task_metrics <<- wr_result %>%
            mutate(task = task_name,
                   dataset_id = dataset_id,
                   ident = ident_col)
        }
      },
      "biological_perturbations" = {
        specs_file <- file.path(proj_root(), "data_processing", "config", "specs_datasets.csv")
        wr_result <<- run_task_biological_perturbations(obj_prepared, config, config[["task_biological_perturbations"]],
                                                        task_output_dir, specs_path = specs_file,
                                                        results_root  = config$output$dir,
                                                        dataset_stems = dataset_stems)
        if (!is.null(wr_result)) {
          task_metrics <<- wr_result %>%
            mutate(task = task_name,
                   dataset_id = dataset_id,
                   ident = ident_col)
        }
      },
      {
        stop("Unknown task: ", task_name)
      }
    )
    
  }, error = function(e) {
    message_step("ERROR", sprintf("Task execution failed: %s", e$message))
  })
  
  # ========== STEP 4: Save results ==========
  if (!is.null(task_metrics)) {
    save_ok <- tryCatch({
      save_task_results(
        results = task_metrics,
        wr_object = wr_result,
        task_name = task_name,
        dataset_id = dataset_id,
        ident = ident_col,
        output_dir = task_output_dir,
        config = config,
        save_wr = config$output$save_wr_objects
      )
      TRUE
    }, error = function(e) {
      message_step("ERROR", sprintf("Saving results failed: %s", e$message))
      FALSE
    })

    if (!save_ok) {
      return(data.frame(
        dataset_id = dataset_id,
        ident = ident_col,
        task = task_name,
        status = "failed",
        error = "save_failed",
        n_results = 0
      ))
    }
    
    # Clean up
    rm(obj, obj_prepared, wr_result)
    gc()
    
    return(task_metrics %>%
      mutate(
        dataset_id = dataset_id,
        ident = ident_col,
        status = "success",
        error = NA,
        .before = 1
      ))
  } else {
    return(data.frame(
      dataset_id = dataset_id,
      ident = ident_col,
      task = task_name,
      status = "failed",
      error = "task_execution_failed",
      n_results = 0
    ))
  }
}

