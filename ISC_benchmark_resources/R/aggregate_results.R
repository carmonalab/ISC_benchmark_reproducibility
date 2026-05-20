#!/usr/bin/env Rscript
# Aggregate resource benchmarking results from all datasets/idents
#
# Usage:
#   Rscript R/aggregate_results.R [--output OUTPUT_FILE]

# Setup
project_root <- Sys.getenv("PROJECT_ROOT", unset = getwd())
setwd(project_root)

source("R/cli_utils.R")
source("R/shared_helpers.R")
source(file.path("ISC_benchmark_resources", "R", "resources_utils.R"))

# Parse arguments
args <- commandArgs(trailingOnly = TRUE)
output_file <- "ISC_benchmark_resources/results/aggregated_benchmarks.rds"

if (length(args) > 0) {
  for (i in seq_along(args)) {
    if (args[i] == "--output" && i < length(args)) {
      output_file <- args[i + 1]
    }
  }
}

# Ensure output directory exists
output_dir <- dirname(output_file)
resource_ensure_dir(output_dir)

# Load configuration
config <- load_resource_config()
output_root <- config$output_root

resource_message_time("Aggregating benchmarking results from: ", output_root)

# Recursively find all .rds files in output structure
results_files <- list.files(
  output_root,
  pattern = "\\.rds$",
  recursive = TRUE,
  full.names = TRUE
)

if (length(results_files) == 0) {
  resource_message_time("WARNING: No result files found in: ", output_root)
  aggregated <- data.frame(
    duration_ms = numeric(),
    peak_memory_MB = numeric(),
    duration = numeric(),
    memory_usage_MB = numeric(),
    cpu_usage = numeric(),
    method = character(),
    consistency_metric = character(),
    dissimilarity_method = character(),
    dataset = character(),
    dataset_id = character(),
    ident = character(),
    nfeatures = integer(),
    ncells = integer(),
    nsamples = integer(),
    sparsity = numeric(),
    gene_list = character(),
    benchmark_ncores = integer(),
    benchmark_iterations = integer(),
    stringsAsFactors = FALSE
  )
  saveRDS(aggregated, output_file)
  resource_message_time("Saved empty aggregated results to: ", output_file)
  quit(status = 0)
}

resource_message_time("Found ", length(results_files), " result files to aggregate")

all_results <- lapply(results_files, function(file) {
  resource_message_time("  Loading: ", basename(file))
  tryCatch({
    readRDS(file)
  }, error = function(e) {
    resource_message_time("    ERROR reading file: ", conditionMessage(e))
    NULL
  })
})

all_results <- Filter(Negate(is.null), all_results)

if (length(all_results) == 0) {
  resource_message_time("ERROR: Could not read any result files")
  quit(status = 1)
}

aggregated <- dplyr::bind_rows(all_results)

duration_col <- if ("duration_ms" %in% names(aggregated)) "duration_ms" else "duration"
memory_col <- if ("peak_memory_MB" %in% names(aggregated)) "peak_memory_MB" else "memory_usage_MB"

resource_message_time(
  "Aggregated ", nrow(aggregated), " benchmark results from ",
  length(unique(aggregated$dataset_id)), " datasets"
)

message("\n=== Summary Statistics ===")
summary_by_dataset <- aggregated %>%
  dplyr::group_by(dataset_id) %>%
  dplyr::summarise(
    n_results = dplyr::n(),
    n_idents = dplyr::n_distinct(ident),
    n_metric_combos = dplyr::n_distinct(
      paste0(dissimilarity_method, "::", consistency_metric)
    ),
    mean_duration_ms = mean(.data[[duration_col]], na.rm = TRUE),
    mean_peak_memory_MB = mean(.data[[memory_col]], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(dataset_id)

print(summary_by_dataset)

saveRDS(aggregated, output_file)
resource_message_time("Saved aggregated results to: ", output_file)

csv_file <- sub("\\.rds$", ".csv", output_file)
write.csv(aggregated, csv_file, row.names = FALSE)
resource_message_time("Saved CSV summary to: ", csv_file)

resource_message_time("Aggregation completed successfully!")
