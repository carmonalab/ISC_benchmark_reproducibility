# label_transfer_task/_targets.R --- Targets workflow for label-transfer benchmark

suppressPackageStartupMessages({
  library(targets)
  library(tarchetypes)
  library(tidyverse)
})

# Load all pipeline-specific functions
source("R/00_utils.R")
source("R/01_classifiers.R")
source("R/02_scTypeEval_helpers.R")
source("R/03_plots_tables.R")

# ============================================================================
# CONFIGURATION LAYER
# ============================================================================

list(
  # Load all parameters from pipeline-specific config
  # (label_transfer_task/config/label_transfer_parameters.yaml contains all required settings)
  tar_target(
    params,
    load_pipeline_config("label_transfer_parameters.yaml"),
    format = "file"
  ),
  
  # Random seed
  tar_target(
    seed_global,
    params$seed %||% 42
  ),
  
  # ========================================================================
  # DATA LAYER: Identify datasets with label_transfer data
  # ========================================================================
  
  tar_target(
    lt_dataset_ids,
    {
      # Find all datasets that have prepared label_transfer data
      data_dir <- lt_data_processed_dir()
      
      if (!dir.exists(data_dir)) {
        warning("Label-transfer data directory not found: ", data_dir)
        return(character(0))
      }
      
      # List all dataset subdirectories that have query.rds + reference.rds
      dataset_dirs <- list.dirs(data_dir, recursive = FALSE, full.names = FALSE)
      
      # Filter by existence of both files
      valid_datasets <- dataset_dirs[
        map_lgl(dataset_dirs, ~ {
          validate_label_transfer_data(.x)
        })
      ]
      
      # Also filter by participation flag (optional)
      valid_datasets[
        map_lgl(valid_datasets, validate_label_transfer_participation)
      ]
    }
  ),
  
  tar_target(
    lt_classifiers,
    get_lt_classifiers(params)
  ),
  
  # ========================================================================
  # GRID EXPANSION: datasets × classifiers × replicates
  # ========================================================================
  
  tar_target(
    lt_grid,
    {
      expand_grid(
        dataset_id = lt_dataset_ids,
        classifier = lt_classifiers,
        replicate = 1:3  # 3 replicates per classifier-dataset pair
      )
    }
  ),
  
  # ========================================================================
  # LABEL-TRANSFER EXECUTION
  # ========================================================================
  
  tar_target(
    lt_classifier_results,
    {
      run_label_transfer_classifier(
        dataset_id = lt_grid$dataset_id,
        classifier_name = lt_grid$classifier,
        rep = lt_grid$replicate,
        data_dir = lt_data_processed_dir(),
        output_dir = lt_raw_results_dir(),
        seed = seed_global + lt_grid$replicate
      )
    },
    format = "file",
    pattern = map(lt_grid),
    iteration = "list"
  ),
  
  # ========================================================================
  # AGGREGATION & RESULTS
  # ========================================================================
  
  tar_target(
    lt_aggregated_results,
    {
      aggregate_label_transfer_results(
        results_dir = lt_raw_results_dir(),
        output_file = file.path(
          lt_aggregated_dir(),
          "label_transfer_metrics_aggregated.csv"
        )
      )
    },
    format = "file"
  ),
  
  # ========================================================================
  # SUMMARY & PLOTTING
  # ========================================================================
  
  tar_target(
    lt_summary_stats,
    {
      results <- read.csv(lt_aggregated_results)
      summary <- summarize_label_transfer_results(results)
      output_file <- file.path(
        lt_aggregated_dir(),
        "label_transfer_summary_stats.csv"
      )
      write.csv(summary, output_file, row.names = FALSE)
      invisible(output_file)
    },
    format = "file"
  ),
  
  tar_target(
    lt_figures,
    {
      results <- read.csv(lt_aggregated_results)
      plot_label_transfer_benchmarks(
        results_table = results,
        output_dir = lt_figures_dir()
      )
    },
    format = "file"
  )
)
