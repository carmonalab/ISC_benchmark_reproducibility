# label_transfer_task/R/03_plots_tables.R --- Aggregation and plotting

source("R/00_utils.R")

# ============================================================================
# AGGREGATE RESULTS
# ============================================================================

aggregate_label_transfer_results <- function(results_dir, output_file = NULL) {
  
  message_step("aggregate_label_transfer_results", 
               "Reading all classifier results...")
  
  # Find all .rds files matching pattern
  result_files <- list.files(
    results_dir,
    pattern = "^.+_rep[0-9]+\\.rds$",
    full.names = TRUE
  )
  
  if (length(result_files) == 0) {
    warning("No label-transfer result files found in ", results_dir)
    return(NULL)
  }
  
  # Read and combine all results
  results_list <- map_df(result_files, ~ {
    df <- readRDS(.x)
    # Add metadata from sidecar if available
    yaml_path <- paste0(.x, ".yaml")
    if (file.exists(yaml_path)) {
      meta <- load_config(yaml_path)
      df$seed <- meta$seed
      df$timestamp <- as.character(meta$timestamp)
    }
    df
  })
  
  message_step("aggregated", nrow(results_list), " total results")
  
  # Save if output path specified
  if (!is.null(output_file)) {
    ensure_dir(dirname(output_file))
    write.csv(results_list, output_file, row.names = FALSE)
    message_step("saved", output_file)
    invisible(output_file)
  } else {
    invisible(results_list)
  }
}

# ============================================================================
# SUMMARY STATISTICS
# ============================================================================

summarize_label_transfer_results <- function(results_table) {
  
  results_table %>%
    group_by(classifier) %>%
    summarise(
      n_datasets = n_distinct(dataset_id),
      n_runs = n(),
      mean_accuracy = mean(accuracy, na.rm = TRUE),
      sd_accuracy = sd(accuracy, na.rm = TRUE),
      mean_balanced_accuracy = mean(balanced_accuracy, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(mean_accuracy))
}

# ============================================================================
# PLACEHOLDER: PLOTTING FUNCTIONS
# ============================================================================

plot_label_transfer_benchmarks <- function(results_table, output_dir) {
  
  ensure_dir(output_dir)
  message_step("plot_label_transfer_benchmarks", 
               "Generating plots...")
  
  # Placeholder: will add ggplot2-based figures
  # Examples:
  #  - Figure 5A: Classifier accuracy by dataset
  #  - Figure 5B: Balanced accuracy comparison
  #  - Figure 5C: per-cell-type confusion matrices
  
  # For now, save a summary
  summary <- summarize_label_transfer_results(results_table)
  write.csv(
    summary,
    file.path(output_dir, "label_transfer_summary.csv"),
    row.names = FALSE
  )
  
  invisible(list.files(output_dir, full.names = TRUE))
}
