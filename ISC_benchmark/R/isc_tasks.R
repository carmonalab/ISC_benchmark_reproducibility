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
  ext_cfg <- get_external_methods_config(config)

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

  build_baseline_external_state <- function(task_name, obj_prepared) {
    base_state <- list(
      rep = 1L,
      original_ident = ident_col,
      perturbed_ctype = NA_character_,
      active_ident = ident_col
    )

    switch(task_name,
      "missclassify" = {
        c(base_state, list(rate = 1))
      },
      "SplitCelltype" = {
        c(base_state, list(rate = 1))
      },
      "Nsamples" = {
        c(base_state, list(rate = 1))
      },
      "NCell" = {
        c(base_state, list(rate = 1))
      },
      "Nct" = {
        all_cts <- unique(obj_prepared$metadata[[obj_prepared$ident]])
        all_cts <- all_cts[!is.na(all_cts)]
        c(base_state, list(rate = paste(sort(as.character(all_cts)), collapse = "-"), rep = NA_integer_))
      },
      "cellular_complexity" = {
        n_cts <- length(unique(obj_prepared$metadata[[obj_prepared$ident]]))
        c(base_state, list(rate = as.numeric(n_cts), rep = NA_integer_))
      },
      NULL
    )
  }

  expected_external_metrics <- function(ext_cfg, task_name) {
    metrics <- character(0)
    if (!is.null(ext_cfg) && is_external_method_enabled_for_task(ext_cfg$sccaf, task_name)) {
      metrics <- c(metrics, "SCCAF")
    }
    if (!is.null(ext_cfg) && is_external_method_enabled_for_task(ext_cfg$anticor_features, task_name)) {
      metrics <- c(metrics, "anticor_features")
    }
    if (!is.null(ext_cfg) && is_external_method_enabled_for_task(ext_cfg$scshc, task_name)) {
      metrics <- c(metrics, "scSHC")
    }
    unique(metrics)
  }

  baseline_external_rows_for_task <- function(baseline_df, task_name, ext_cfg, obj_prepared) {
    if (is.null(baseline_df) || is.null(ext_cfg) || !(task_name %in% TASKS_WITH_BASELINE)) {
      return(NULL)
    }

    needed_metrics <- expected_external_metrics(ext_cfg, task_name)
    if (length(needed_metrics) == 0 || !"consistency_metric" %in% names(baseline_df)) {
      return(NULL)
    }

    external_bl <- baseline_df %>%
      dplyr::filter(consistency_metric %in% needed_metrics)

    if (nrow(external_bl) == 0) {
      return(NULL)
    }

    switch(task_name,
      "missclassify" = baseline_for_task(external_bl, "missclassify", filter_external = FALSE),
      "SplitCelltype" = baseline_for_task(external_bl, "SplitCelltype", filter_external = FALSE),
      "Nsamples" = baseline_for_task(external_bl, "Nsamples", filter_external = FALSE),
      "NCell" = baseline_for_task(external_bl, "NCell", filter_external = FALSE),
      "Nct" = {
        all_cts <- unique(obj_prepared$metadata[[obj_prepared$ident]])
        all_cts <- all_cts[!is.na(all_cts)]
        baseline_for_Nct(external_bl, all_cts, filter_external = FALSE)
      },
      "cellular_complexity" = {
        n_cts <- length(unique(obj_prepared$metadata[[obj_prepared$ident]]))
        baseline_for_mergeCT(external_bl, n_cts, filter_external = FALSE)
      },
      NULL
    )
  }

  if (!is.null(baseline_df) && !is.null(ext_cfg) && task_name %in% TASKS_WITH_BASELINE) {
    needed_metrics <- expected_external_metrics(ext_cfg, task_name)
    existing_metrics <- if ("consistency_metric" %in% names(baseline_df)) {
      unique(as.character(stats::na.omit(baseline_df$consistency_metric)))
    } else {
      character(0)
    }
    missing_metrics <- setdiff(needed_metrics, existing_metrics)

    if (length(missing_metrics) > 0) {
      baseline_state <- build_baseline_external_state(task_name, obj_prepared)
      baseline_external <- tryCatch(
        {
          sc_baseline_raw <- scTypeEval::create_scTypeEval(
            matrix = obj_prepared$count_matrix,
            metadata = obj_prepared$metadata,
            active_ident = ident_col
          )

          sc_baseline <- scTypeEval::wrapper_scTypeEval(
            sc_baseline_raw,
            ident = ident_col,
            sample = config$common$sample,
            gene_list = NULL,
            reduction = config$common$reduction,
            ndim = config$common$ndim,
            normalization_method = config$common$normalization_method,
            dissimilarity_method = config$common$dissimilarity_method,
            min_samples = config$common$min_samples,
            min_cells = config$common$min_cells,
            verbose = isTRUE(config$common$verbose)
          )

          run_external_methods_for_state(
            sc_obj = sc_baseline,
            task_name = task_name,
            dataset_id = dataset_id,
            ident_col = ident_col,
            output_dir = file.path(output_dir, sprintf("%s_%s", task_name, ident_col)),
            config = config,
            state = baseline_state
          )
        },
        error = function(e) {
          message_step("EXTERNAL", sprintf("Baseline external cache enrichment failed for %s: %s", task_name, e$message))
          NULL
        }
      )

      if (!is.null(baseline_external) && nrow(baseline_external) > 0) {
        baseline_external$task <- "Baseline"
        baseline_df <- dplyr::bind_rows(baseline_df, baseline_external) %>% dplyr::distinct()
        baseline_cache_path <- file.path(output_dir, paste0("baseline_isc_", ident_col, ".rds"))
        saveRDS(baseline_df, baseline_cache_path)
        message_step("BASELINE", sprintf("Cached baseline external rows to %s", basename(baseline_cache_path)))
      }
    }
  }

  # ========== STEP 3: Execute task ==========
  task_output_dir <- file.path(output_dir, sprintf("%s_%s", task_name, ident_col))
  metrics_file <- file.path(
    task_output_dir,
    sprintf("%s_%s_%s_metrics.rds", dataset_id, task_name, ident_col)
  )

  if (file.exists(metrics_file)) {
    message_step("TASK", sprintf("Resuming from saved results: %s", basename(metrics_file)))
    cached_metrics <- readRDS(metrics_file)
    if (!all(c("dataset_id", "status", "error") %in% names(cached_metrics))) {
      cached_metrics <- cached_metrics %>%
        mutate(
          dataset_id = dataset_id,
          ident = ident_col,
          status = "success",
          error = NA,
          .before = 1
        )
    }
    return(cached_metrics)
  }
  
  wr_result <- NULL
  task_metrics <- NULL
  skip_persist <- FALSE
  external_state_callback <- NULL

  if (!is.null(ext_cfg)) {
    external_state_callback <- function(sc_obj_state, state) {
      run_external_methods_for_state(
        sc_obj = sc_obj_state,
        task_name = task_name,
        dataset_id = dataset_id,
        ident_col = ident_col,
        output_dir = task_output_dir,
        config = config,
        state = state
      )
    }
  }
  
  tryCatch({
    message_step("TASK", sprintf("Running %s...", task_name))
    
    # Extract task-specific configuration from main config
    task_config_key <- paste0("task_", task_name)
    task_config <- config[[task_config_key]][["isc_params"]]
    
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
        wr_result <- run_task_missclassify(obj_prepared, config, task_config, task_output_dir,
                                             baseline_df = baseline_df,
                                             external_state_callback = external_state_callback)
        task_metrics <- wr_result
      },
      "SplitCelltype" = {
        wr_result <- run_task_SplitCelltype(obj_prepared, config, task_config, task_output_dir,
                                              baseline_df = baseline_df,
                                              external_state_callback = external_state_callback)
        task_metrics <- wr_result
      },
      "Nct" = {
        wr_result <- run_task_Nct(obj_prepared, config, task_config, task_output_dir,
                                    baseline_df = baseline_df,
                                    external_state_callback = external_state_callback)
        task_metrics <- wr_result
      },
      "cellular_complexity" = {
        wr_result <- run_task_cellular_complexity(obj_prepared, config, task_config, task_output_dir,
                                                   baseline_df = baseline_df,
                                                   external_state_callback = external_state_callback)
        task_metrics <- wr_result
      },
      "Nsamples" = {
        wr_result <- run_task_Nsamples(obj_prepared, config, task_config, task_output_dir,
                                         baseline_df = baseline_df,
                                         external_state_callback = external_state_callback)
        task_metrics <- wr_result
      },
      "NCell" = {
        wr_result <- run_task_NCell(obj_prepared, config, task_config, task_output_dir,
                                      baseline_df = baseline_df,
                                      external_state_callback = external_state_callback)
        task_metrics <- wr_result
      },
      "batch_effects" = {
        specs_file <- file.path(proj_root(), "data_processing", "config", "specs_datasets.csv")
        wr_result <- run_task_batch_effects(obj_prepared, config, config[["task_batch_effects"]], task_output_dir,
                                             specs_path    = specs_file,
                                             results_root  = config$output$dir,
                                             dataset_stems = dataset_stems,
                                             external_state_callback = external_state_callback)
        if (!is.null(wr_result)) {
          task_metrics <- wr_result %>%
            mutate(task = task_name,
                   dataset_id = dataset_id,
                   ident = ident_col)
        } else {
          # No resolvable batch pairs is an expected outcome for some families.
          # Mark as success but skip persistence so no files/directories are created.
          skip_persist <- TRUE
        }
      },
      "biological_perturbations" = {
        specs_file <- file.path(proj_root(), "data_processing", "config", "specs_datasets.csv")
        wr_result <- run_task_biological_perturbations(obj_prepared, config, config[["task_biological_perturbations"]],
                                                        task_output_dir, specs_path = specs_file,
                                                        results_root  = config$output$dir,
                                                        dataset_stems = dataset_stems,
                                                        external_state_callback = external_state_callback)
        if (!is.null(wr_result)) {
          task_metrics <- wr_result %>%
            mutate(task = task_name,
                   dataset_id = dataset_id,
                   ident = ident_col)
        } else {
          # No resolvable perturbation pairs is an expected outcome for some families.
          # Mark as success but skip persistence so no files/directories are created.
          skip_persist <- TRUE
        }
      },
      {
        stop("Unknown task: ", task_name)
      }
    )
    
  }, error = function(e) {
    message_step("ERROR", sprintf("Task execution failed: %s", e$message))
  })

  if (isTRUE(skip_persist)) {
    return(data.frame(
      dataset_id = dataset_id,
      ident = ident_col,
      task = task_name,
      status = "success",
      error = NA,
      n_results = 0,
      stringsAsFactors = FALSE
    ))
  }
  
  # ========== STEP 4: Save results ==========
  if (!is.null(task_metrics)) {
    external_metrics <- attr(wr_result, "external_state_scores")
    task_metrics <- coerce_rep_to_character(task_metrics)
    external_metrics <- coerce_rep_to_character(external_metrics)

    if (!is.null(ext_cfg) && !is.null(baseline_df) && task_name %in% TASKS_WITH_BASELINE) {
      baseline_external_metrics <- baseline_external_rows_for_task(baseline_df, task_name, ext_cfg, obj_prepared)
      baseline_external_metrics <- coerce_rep_to_character(baseline_external_metrics)
      if (!is.null(baseline_external_metrics) && nrow(baseline_external_metrics) > 0) {
        external_metrics <- dplyr::bind_rows(baseline_external_metrics, external_metrics)
      }
    }

    if (is.null(external_metrics) && !is.null(ext_cfg)) {
      external_metrics <- tryCatch(
        run_external_methods_for_task(
          obj_prepared = obj_prepared,
          task_metrics = task_metrics,
          task_name = task_name,
          dataset_id = dataset_id,
          ident_col = ident_col,
          output_dir = task_output_dir,
          config = config
        ),
        error = function(e) {
          message_step("EXTERNAL", sprintf("External methods failed for %s: %s", task_name, e$message))
          NULL
        }
      )
    }

    if (!is.null(external_metrics) && nrow(external_metrics) > 0) {
      task_metrics <- dplyr::bind_rows(task_metrics, external_metrics)
    }

    if (!("level" %in% colnames(task_metrics))) {
      task_metrics$level <- "celltype"
    } else {
      task_metrics$level <- ifelse(
        !is.na(task_metrics$celltype) & task_metrics$celltype == "__global__",
        "global",
        ifelse(is.na(task_metrics$level) | task_metrics$level == "", "celltype", as.character(task_metrics$level))
      )
    }

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

