# label_transfer_task/_targets.R --- Targets workflow for label-transfer benchmark

suppressPackageStartupMessages({
  library(targets)
  library(tarchetypes)
  library(tidyverse)
})

# Shared messaging helpers (message_time/message_step)
source("../R/shared_helpers.R")
source(proj_path("R/cli_utils.R"))

# Load all pipeline-specific functions
source("R/00_utils.R")
source("R/01_classifiers.R")
source("R/02_plots_tables.R")
source("R/03_consistency.R")

# ============================================================================
# CONFIGURATION LAYER
# ============================================================================

list(
  # Load all parameters from pipeline-specific config
  # (label_transfer_task/config/label_transfer_parameters.yaml contains all required settings)
  tar_target(
    lt_params,
    load_pipeline_config("label_transfer_parameters.yaml")
  ),
  
  # Random seed
  tar_target(
    seed_global,
    get_lt_seed(lt_params)
  ),

  tar_target(
    n_cores,
    get_lt_n_cores(lt_params)
  ),

  tar_target(
    n_replicates,
    get_lt_n_replicates(lt_params)
  ),
  
  # ========================================================================
  # DATA LAYER: Identify datasets with label_transfer data
  # ========================================================================
  
  tar_target(
    lt_dataset_ids,
    {
      list_label_transfer_datasets_from_isc(lt_params)
    }
  ),
  
  tar_target(
    lt_classifiers,
    get_lt_classifiers(lt_params)
  ),
  
  # ========================================================================
  # GRID EXPANSION: datasets × classifiers × replicates
  # ========================================================================

  tar_target(
    lt_split_grid,
    {
      expand_grid(
        dataset_id = lt_dataset_ids,
        replicate = seq_len(n_replicates)
      )
    }
  ),

  tar_target(
    lt_prepared_splits,
    {
      prepare_label_transfer_split(
        dataset_id = lt_split_grid$dataset_id,
        replicate = lt_split_grid$replicate,
        params = lt_params
      )
    },
    format = "file",
    pattern = map(lt_split_grid),
    iteration = "list"
  ),
  
  tar_target(
    lt_grid,
    {
      expand_grid(
        dataset_id = lt_dataset_ids,
        classifier = lt_classifiers,
        replicate = seq_len(n_replicates)
      )
    }
  ),
  
  # ========================================================================
  # LABEL-TRANSFER EXECUTION
  # ========================================================================
  
  tar_target(
    lt_classifier_results,
    {
      invisible(lt_prepared_splits)
      run_label_transfer_classifier(
        dataset_id = lt_grid$dataset_id,
        classifier_name = lt_grid$classifier,
        rep = lt_grid$replicate,
        data_dir = lt_data_processed_dir(lt_grid$replicate),
        output_dir = lt_raw_results_dir(),
        seed = seed_global + lt_grid$replicate - 1,
        ncores = n_cores
      )
    },
    format = "file",
    pattern = map(lt_grid),
    iteration = "list"
  ),

  # ========================================================================
  # CONSISTENCY METRICS
  #   - Query split: consistency on predicted labels + F1 vs ground truth
  #   - Reference split: consistency on ground truth only
  # ========================================================================

  tar_target(
    lt_query_consistency,
    {
      compute_lt_query_consistency(
        dataset_id = lt_grid$dataset_id,
        classifier_name = lt_grid$classifier,
        rep = lt_grid$replicate,
        result_path = lt_classifier_results,
        data_dir = lt_data_processed_dir(lt_grid$replicate),
        output_dir = lt_consistency_dir(),
        ncores = n_cores
      )
    },
    format = "file",
    pattern = map(lt_grid, lt_classifier_results),
    iteration = "list"
  ),

  tar_target(
    lt_reference_consistency,
    {
      compute_lt_reference_consistency(
        dataset_id = lt_split_grid$dataset_id,
        rep = lt_split_grid$replicate,
        data_dir = lt_data_processed_dir(lt_split_grid$replicate),
        output_dir = lt_consistency_dir(),
        ncores = n_cores
      )
    },
    format = "file",
    pattern = map(lt_split_grid),
    iteration = "list"
  ),

  tar_target(
    lt_consistency_aggregated,
    {
      list(lt_query_consistency, lt_reference_consistency)
      out <- aggregate_lt_consistency_results(
        consistency_dir = lt_consistency_dir(),
        output_file = file.path(
          lt_aggregated_dir(),
          "label_transfer_consistency.csv"
        )
      )
      if (is.null(out)) character(0) else out
    },
    format = "file"
  ),
  
  # ========================================================================
  # AGGREGATION & RESULTS
  # ========================================================================
  
  tar_target(
    lt_aggregated_results,
    {
      list(lt_classifier_results)
      out <- aggregate_label_transfer_results(
        results_dir = lt_raw_results_dir(),
        output_file = file.path(
          lt_aggregated_dir(),
          "label_transfer_metrics_aggregated.csv"
        )
      )
      if (is.null(out)) character(0) else out
    },
    format = "file"
  ),
  
  # ========================================================================
  # SUMMARY & PLOTTING
  # ========================================================================
  
  tar_target(
    lt_summary_stats,
    {
      if (is.null(lt_aggregated_results) || !file.exists(lt_aggregated_results)) {
        message("[SKIP] lt_summary_stats: no aggregated results available")
        return(character(0))
      }
      results <- read.csv(lt_aggregated_results)
      summary <- summarize_label_transfer_results(results)
      output_file <- file.path(
        lt_aggregated_dir(),
        "label_transfer_summary_stats.csv"
      )
      write.csv(summary, output_file, row.names = FALSE)
      output_file
    },
    format = "file"
  ),
  
  tar_target(
    lt_figures,
    {
      if (is.null(lt_aggregated_results) || !file.exists(lt_aggregated_results)) {
        message("[SKIP] lt_figures: no aggregated results available")
        return(character(0))
      }
      results <- read.csv(lt_aggregated_results)
      plot_label_transfer_benchmarks(
        results_table = results,
        output_dir = lt_figures_dir()
      )
    },
    format = "file"
  )
)
