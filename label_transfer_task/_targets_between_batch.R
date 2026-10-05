# label_transfer_task/_targets_between_batch.R --- Targets workflow for batch label transfer
# (one reference dataset vs. all other same-annotation datasets combined as query)

suppressPackageStartupMessages({
  library(targets)
  library(tarchetypes)
  library(tidyverse)
})

source("../R/shared_helpers.R")
source(proj_path("R/cli_utils.R"))

source("R/00_utils.R")
source("R/00_between_utils.R")
source("R/00_batch_utils.R")
source("R/01_classifiers.R")
source("R/02_plots_tables.R")
source("R/03_consistency.R")

tar_option_set(error = "continue")

list(
  tar_target(
    lt_batch_params,
    load_pipeline_config("label_transfer_batch_parameters.yaml")
  ),

  tar_target(
    lt_batch_seed,
    get_lt_seed(lt_batch_params)
  ),

  tar_target(
    lt_batch_n_cores,
    get_lt_n_cores(lt_batch_params)
  ),

  tar_target(
    lt_batch_n_replicates,
    get_lt_n_replicates(lt_batch_params)
  ),

  tar_target(
    lt_batch_classifiers,
    get_lt_classifiers(lt_batch_params)
  ),

  tar_target(
    lt_batch_pairs,
    list_batch_pairs(lt_batch_params)
  ),

  tar_target(
    lt_batch_prepared_pairs,
    prepare_batch_pair(
      pair_id = lt_batch_pairs$pair_id,
      reference_dataset_id = lt_batch_pairs$reference_dataset_id,
      query_dataset_ids = lt_batch_pairs$query_dataset_ids,
      query_dataset_id = lt_batch_pairs$query_dataset_id
    ),
    format = "file",
    pattern = map(lt_batch_pairs),
    iteration = "list"
  ),

  tar_target(
    lt_batch_unique_ref_grid,
    tidyr::crossing(
      reference_dataset_id = unique(lt_batch_pairs$reference_dataset_id),
      replicate = seq_len(lt_batch_n_replicates)
    )
  ),

  tar_target(
    lt_batch_grid,
    tidyr::crossing(
      lt_batch_pairs,
      classifier = lt_batch_classifiers,
      replicate = seq_len(lt_batch_n_replicates)
    )
  ),

  tar_target(
    lt_batch_ensemble_grid,
    tidyr::crossing(
      lt_batch_pairs,
      replicate = seq_len(lt_batch_n_replicates)
    )
  ),

  tar_target(
    lt_batch_classifier_results,
    {
      invisible(lt_batch_prepared_pairs)
      run_label_transfer_classifier(
        dataset_id = lt_batch_grid$pair_id,
        classifier_name = lt_batch_grid$classifier,
        rep = lt_batch_grid$replicate,
        data_dir = lt_batch_data_processed_dir(),
        output_dir = lt_batch_raw_results_dir(),
        seed = lt_batch_seed + lt_batch_grid$replicate - 1,
        ncores = lt_batch_n_cores,
        reference_dataset_id = lt_batch_grid$reference_dataset_id,
        query_dataset_id = lt_batch_grid$query_dataset_id
      )
    },
    format = "file",
    pattern = map(lt_batch_grid),
    iteration = "list"
  ),

  tar_target(
    lt_batch_ensemble_results,
    {
      invisible(lt_batch_classifier_results)
      
      run_ensemble_classifier_between_datasets(
        pair_id = lt_batch_ensemble_grid$pair_id,
        reference_dataset_id = lt_batch_ensemble_grid$reference_dataset_id,
        query_dataset_id = lt_batch_ensemble_grid$query_dataset_id,
        rep = lt_batch_ensemble_grid$replicate,
        output_dir = lt_batch_raw_results_dir()
      )
    },
    format = "file",
    pattern = map(lt_batch_ensemble_grid),
    iteration = "list"
  ),

  tar_target(
    lt_batch_query_consistency,
    {
      compute_lt_query_consistency(
        dataset_id = lt_batch_grid$pair_id,
        classifier_name = lt_batch_grid$classifier,
        rep = lt_batch_grid$replicate,
        result_path = lt_batch_classifier_results,
        data_dir = lt_batch_data_processed_dir(),
        output_dir = lt_batch_consistency_dir(),
        ncores = lt_batch_n_cores,
        reference_dataset_id = lt_batch_grid$reference_dataset_id,
        query_dataset_id = lt_batch_grid$query_dataset_id
      )
    },
    format = "file",
    pattern = map(lt_batch_grid, lt_batch_classifier_results),
    iteration = "list"
  ),

  tar_target(
    lt_batch_ensemble_consistency,
    {
      compute_lt_query_consistency(
        dataset_id = lt_batch_ensemble_grid$pair_id,
        classifier_name = "Ensemble",
        rep = lt_batch_ensemble_grid$replicate,
        result_path = lt_batch_ensemble_results,
        data_dir = lt_batch_data_processed_dir(),
        output_dir = lt_batch_consistency_dir(),
        ncores = lt_batch_n_cores,
        reference_dataset_id = lt_batch_ensemble_grid$reference_dataset_id,
        query_dataset_id = lt_batch_ensemble_grid$query_dataset_id
      )
    },
    format = "file",
    pattern = map(lt_batch_ensemble_grid, lt_batch_ensemble_results),
    iteration = "list"
  ),

  tar_target(
    lt_batch_unique_ref_cons,
    {
      invisible(lt_batch_prepared_pairs)
      compute_lt_batch_unique_reference_consistency(
        reference_dataset_id = lt_batch_unique_ref_grid$reference_dataset_id,
        rep = lt_batch_unique_ref_grid$replicate,
        pairs = lt_batch_pairs,
        data_dir = lt_batch_data_processed_dir(),
        ncores = lt_batch_n_cores
      )
    },
    pattern = map(lt_batch_unique_ref_grid),
    iteration = "list"
  ),

  tar_target(
    lt_batch_reference_consistency,
    {
      write_lt_batch_reference_consistency_outputs(
        unique_cons = lt_batch_unique_ref_cons,
        unique_grid = lt_batch_unique_ref_grid,
        pairs = lt_batch_pairs,
        output_dir = lt_batch_consistency_dir()
      )
    },
    format = "file"
  ),

  tar_target(
    lt_batch_consistency_aggregated,
    {
      list(lt_batch_query_consistency, lt_batch_ensemble_consistency, lt_batch_reference_consistency)
      out <- aggregate_lt_consistency_results(
        consistency_dir = lt_batch_consistency_dir(),
        output_file = file.path(
          lt_batch_aggregated_dir(),
          "label_transfer_batch_consistency.csv"
        )
      )
      if (is.null(out)) character(0) else out
    },
    format = "file"
  ),

  tar_target(
    lt_batch_aggregated_results,
    {
      list(lt_batch_classifier_results)
      out <- aggregate_label_transfer_results(
        results_dir = lt_batch_raw_results_dir(),
        output_file = file.path(
          lt_batch_aggregated_dir(),
          "label_transfer_batch_metrics_aggregated.csv"
        )
      )
      if (is.null(out)) character(0) else out
    },
    format = "file"
  ),

  tar_target(
    lt_batch_summary_stats,
    {
      if (is.null(lt_batch_aggregated_results) || !file.exists(lt_batch_aggregated_results)) {
        message("[SKIP] lt_batch_summary_stats: no aggregated results available")
        return(character(0))
      }
      results <- read.csv(lt_batch_aggregated_results)
      summary <- summarize_label_transfer_results(results)
      output_file <- file.path(
        lt_batch_aggregated_dir(),
        "label_transfer_batch_summary_stats.csv"
      )
      write.csv(summary, output_file, row.names = FALSE)
      output_file
    },
    format = "file"
  ),

  tar_target(
    lt_batch_figures,
    {
      if (is.null(lt_batch_aggregated_results) || !file.exists(lt_batch_aggregated_results)) {
        message("[SKIP] lt_batch_figures: no aggregated results available")
        return(character(0))
      }
      results <- read.csv(lt_batch_aggregated_results)
      plot_label_transfer_benchmarks(
        results_table = results,
        output_dir = lt_batch_figures_dir()
      )
    },
    format = "file"
  )
)
