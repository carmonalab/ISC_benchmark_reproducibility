# label_transfer_task/_targets_between_datasets.R --- Targets workflow for between-dataset label transfer

suppressPackageStartupMessages({
  library(targets)
  library(tarchetypes)
  library(tidyverse)
})

source("../R/shared_helpers.R")
source(proj_path("R/cli_utils.R"))

source("R/00_utils.R")
source("R/00_between_utils.R")
source("R/01_classifiers.R")
source("R/02_plots_tables.R")
source("R/03_consistency.R")

list(
  tar_target(
    lt_between_params,
    load_pipeline_config("label_transfer_between_parameters.yaml")
  ),

  tar_target(
    lt_between_seed,
    get_lt_seed(lt_between_params)
  ),

  tar_target(
    lt_between_n_cores,
    get_lt_n_cores(lt_between_params)
  ),

  tar_target(
    lt_between_n_replicates,
    get_lt_n_replicates(lt_between_params)
  ),

  tar_target(
    lt_between_classifiers,
    get_lt_classifiers(lt_between_params)
  ),

  tar_target(
    lt_between_pairs,
    list_between_dataset_pairs(lt_between_params)
  ),

  tar_target(
    lt_between_prepared_pairs,
    prepare_between_dataset_pair(
      pair_id = lt_between_pairs$pair_id,
      reference_dataset_id = lt_between_pairs$reference_dataset_id,
      query_dataset_id = lt_between_pairs$query_dataset_id,
      params = lt_between_params
    ),
    format = "file",
    pattern = map(lt_between_pairs),
    iteration = "list"
  ),

  tar_target(
    lt_between_ref_grid,
    tidyr::crossing(
      lt_between_pairs,
      replicate = seq_len(lt_between_n_replicates)
    )
  ),

  tar_target(
    lt_between_grid,
    tidyr::crossing(
      lt_between_pairs,
      classifier = lt_between_classifiers,
      replicate = seq_len(lt_between_n_replicates)
    )
  ),

  tar_target(
    lt_between_classifier_results,
    {
      invisible(lt_between_prepared_pairs)
      run_label_transfer_classifier(
        dataset_id = lt_between_grid$pair_id,
        classifier_name = lt_between_grid$classifier,
        rep = lt_between_grid$replicate,
        data_dir = lt_between_data_processed_dir(),
        output_dir = lt_between_raw_results_dir(),
        seed = lt_between_seed + lt_between_grid$replicate - 1,
        ncores = lt_between_n_cores,
        reference_dataset_id = lt_between_grid$reference_dataset_id,
        query_dataset_id = lt_between_grid$query_dataset_id
      )
    },
    format = "file",
    pattern = map(lt_between_grid),
    iteration = "list"
  ),

  tar_target(
    lt_between_query_consistency,
    {
      compute_lt_query_consistency(
        dataset_id = lt_between_grid$pair_id,
        classifier_name = lt_between_grid$classifier,
        rep = lt_between_grid$replicate,
        result_path = lt_between_classifier_results,
        data_dir = lt_between_data_processed_dir(),
        output_dir = lt_between_consistency_dir(),
        ncores = lt_between_n_cores,
        reference_dataset_id = lt_between_grid$reference_dataset_id,
        query_dataset_id = lt_between_grid$query_dataset_id
      )
    },
    format = "file",
    pattern = map(lt_between_grid, lt_between_classifier_results),
    iteration = "list"
  ),

  tar_target(
    lt_between_reference_consistency,
    {
      compute_lt_reference_consistency(
        dataset_id = lt_between_ref_grid$pair_id,
        rep = lt_between_ref_grid$replicate,
        data_dir = lt_between_data_processed_dir(),
        output_dir = lt_between_consistency_dir(),
        ncores = lt_between_n_cores,
        reference_dataset_id = lt_between_ref_grid$reference_dataset_id,
        query_dataset_id = lt_between_ref_grid$query_dataset_id
      )
    },
    format = "file",
    pattern = map(lt_between_ref_grid),
    iteration = "list"
  ),

  tar_target(
    lt_between_consistency_aggregated,
    {
      list(lt_between_query_consistency, lt_between_reference_consistency)
      out <- aggregate_lt_consistency_results(
        consistency_dir = lt_between_consistency_dir(),
        output_file = file.path(
          lt_between_aggregated_dir(),
          "label_transfer_between_consistency.csv"
        )
      )
      if (is.null(out)) character(0) else out
    },
    format = "file"
  ),

  tar_target(
    lt_between_aggregated_results,
    {
      list(lt_between_classifier_results)
      out <- aggregate_label_transfer_results(
        results_dir = lt_between_raw_results_dir(),
        output_file = file.path(
          lt_between_aggregated_dir(),
          "label_transfer_between_metrics_aggregated.csv"
        )
      )
      if (is.null(out)) character(0) else out
    },
    format = "file"
  ),

  tar_target(
    lt_between_summary_stats,
    {
      if (is.null(lt_between_aggregated_results) || !file.exists(lt_between_aggregated_results)) {
        message("[SKIP] lt_between_summary_stats: no aggregated results available")
        return(character(0))
      }
      results <- read.csv(lt_between_aggregated_results)
      summary <- summarize_label_transfer_results(results)
      output_file <- file.path(
        lt_between_aggregated_dir(),
        "label_transfer_between_summary_stats.csv"
      )
      write.csv(summary, output_file, row.names = FALSE)
      output_file
    },
    format = "file"
  ),

  tar_target(
    lt_between_figures,
    {
      if (is.null(lt_between_aggregated_results) || !file.exists(lt_between_aggregated_results)) {
        message("[SKIP] lt_between_figures: no aggregated results available")
        return(character(0))
      }
      results <- read.csv(lt_between_aggregated_results)
      plot_label_transfer_benchmarks(
        results_table = results,
        output_dir = lt_between_figures_dir()
      )
    },
    format = "file"
  )
)
