#!/usr/bin/env Rscript

#' Prepare Label Transfer Query/Reference Splits
#'
#' Post-processing script to create query (70% samples) and reference (30% samples) 
#' splits from ISC-processed Seurat objects, following:
#' 
#'   Consistency_metrics_benchmark/OPS/manual_classifiers/process_datasets.R
#'
#' Key features:
#'   - Loads Seurat objects from data/processed/isc/
#'   - Extracts counts and metadata using scTypeEval::load_singleCell_object
#'   - Stratified sample-level split (deterministic seed=22)
#'   - Skips datasets with <10 samples
#'   - Outputs query.rds and reference.rds to data/processed/label_transfer/<dataset_id>/
#'
#' Usage:
#'   Rscript label_transfer_task/R/00_prepare_splits.R
#'
#' Requirements:
#'   - Run from project root: cd /path/to/ISC_benchmark_reproducibility
#'   - ISC processed datasets already in data/processed/isc/
#'   - R packages: dplyr, scTypeEval, Seurat, BiocParallel
#'

suppressPackageStartupMessages({
  library(dplyr)
  library(Seurat)
  library(BiocParallel)
  library(yaml)
})

# Source shared utilities (these were previously defined locally here)
source("../../R/cli_utils.R")       # For message_time() and message_step()
source("../../R/shared_helpers.R")  # For ensure_dir(), proj_path(), etc.
source("00_utils.R")                # For get_idents_by_prefix()

# ============================================================================
# Configuration
# ============================================================================

# Fixed seed for reproducible sample splits (matches proc_data.Rmd)
SEED <- 22

# Proportion of samples for query set
PROP_QUERY <- 0.7

# Minimum samples required to process dataset
MIN_SAMPLES <- 10

# ============================================================================
# Helper Functions (module-specific)
# ============================================================================

#' Get identification column for dataset based on filename prefix
#'
#' Maps dataset prefixes to their cell-type identification columns,
#' matching reference implementation exactly.
#' Uses centralized get_idents_by_prefix() function from 00_utils.R
#'
#' @param dataset_key The dataset identifier/filename (e.g., "Stephenson_COVID")
#' @param metadata_cols Available columns in metadata
#'
#' @return Character string: column name to use as cell_type annotation
get_ident_for_dataset <- function(dataset_key, metadata_cols) {
  # Parse dataset_key to infer dataset source
  prefix <- strsplit(dataset_key, "_")[[1]][1]
  
  # Use centralized mapping function
  ident_col <- get_idents_by_prefix(prefix)
  
  # If inferred column not available, try common alternatives
  if (!ident_col %in% metadata_cols) {
    candidates <- c(
      "OriginalAnnotationLevel2", "OriginalAnnotationLevel1",
      "cell_type", "celltype", "cell.type",
      "annotation", "cell_annotation", "type"
    )
    ident_col <- candidates[candidates %in% metadata_cols][1]
    if (is.na(ident_col)) {
      stop("No suitable cell-type column found in metadata. Available: ",
           paste(metadata_cols, collapse = ", "))
    }
  }
  
  ident_col
}
           paste(metadata_cols, collapse = ", "))
    }
  }
  
  ident_col
}

process_label_transfer_dataset <- function(isc_path,
                                          yaml_path,
                                          out_lt_dir) {
  dataset_key <- tools::file_path_sans_ext(basename(isc_path))
  
  message_time(paste("Loading:", dataset_key))
  
  # Load ISC processed object
  obj <- readRDS(isc_path)
  
  if (!inherits(obj, "Seurat")) {
    warning("Not a Seurat object: ", isc_path)
    return(invisible(NULL))
  }
  
  # Check if dataset is marked for label transfer
  if (file.exists(yaml_path)) {
    meta <- yaml::read_yaml(yaml_path)
    if (!isTRUE(meta$label_transfer)) {
      message_time(paste("Skipping (not marked for label transfer):", dataset_key))
      return(invisible(NULL))
    }
  }
  
  # Extract counts and metadata
  counts <- SeuratObject::GetAssayData(obj, assay = DefaultAssay(obj), slot = "counts")
  md <- obj@meta.data
  
  # Check minimum samples threshold
  n_samples <- n_distinct(md$sample)
  if (n_samples < MIN_SAMPLES) {
    message_time(paste("Skipping (n_samples =", n_samples, "<", MIN_SAMPLES, "):", dataset_key))
    return(invisible(NULL))
  }
  
  # Determine cell-type column
  ident_col <- get_ident_for_dataset(dataset_key, colnames(md))
  message_time(paste("Using cell-type column:", ident_col))
  
  # Standardize metadata (matching Consistency_metrics_benchmark/OPS/manual_classifiers)
  md_standardized <- md %>%
    dplyr::select(sample) %>%
    dplyr::mutate(
      cell_type = md[[ident_col]],
      dataset_id = dataset_key
    )
  
  # Split samples: 70% query, 30% reference
  all_samples <- unique(md_standardized$sample)
  n_query <- floor(length(all_samples) * PROP_QUERY)
  
  set.seed(SEED)
  query_samples <- sample(all_samples, size = n_query)
  
  # Subset cells for each split
  cells_query <- rownames(md_standardized)[md_standardized$sample %in% query_samples]
  cells_ref <- rownames(md_standardized)[!md_standardized$sample %in% query_samples]
  
  # Create output directory
  outdir <- ensure_dir(file.path(out_lt_dir, dataset_key))
  
  # Create and save query split
  query <- list(
    counts = counts[, cells_query, drop = FALSE],
    metadata = md_standardized[cells_query, , drop = FALSE]
  )
  saveRDS(query, file.path(outdir, "query.rds"))
  
  # Create and save reference split
  reference <- list(
    counts = counts[, cells_ref, drop = FALSE],
    metadata = md_standardized[cells_ref, , drop = FALSE]
  )
  saveRDS(reference, file.path(outdir, "reference.rds"))
  
  # Summary statistics
  n_query_samples <- n_distinct(query$metadata$sample)
  n_ref_samples <- n_distinct(reference$metadata$sample)
  
  message_time(sprintf(
    "✓ Saved: %s (query: %d samples, %d cells | ref: %d samples, %d cells)",
    outdir,
    n_query_samples, ncol(query$counts),
    n_ref_samples, ncol(reference$counts)
  ))
  
  invisible(list(
    dataset_key = dataset_key,
    outdir = outdir,
    n_query_samples = n_query_samples,
    n_ref_samples = n_ref_samples
  ))
}

main <- function() {
  # Verify we're in correct directory
  if (!file.exists("data_processing/config/datasets.yaml") || !dir.exists("data/processed/isc")) {
    stop("Run from project root: cd /path/to/ISC_benchmark_reproducibility")
  }
  
  message_time("=".repeat(70))
  message_time("Preparing Label Transfer Query/Reference Splits")
  message_time("=".repeat(70))
  
  isc_dir <- "data/processed/isc"
  out_lt_dir <- "data/processed/label_transfer"
  
  ensure_dir(out_lt_dir)
  
  # Find ISC processed RDS files
  isc_files <- list.files(
    isc_dir,
    pattern = "\\.rds$",
    ignore.case = TRUE,
    full.names = TRUE
  )
  
  if (length(isc_files) == 0) {
    stop("No RDS files found in ", isc_dir)
  }
  
  message_time(paste("Processing", length(isc_files), "ISC datasets..."))
  
  # Process all datasets (parallel where possible)
  n_cores <- min(8, length(isc_files))
  
  results <- bplapply(
    isc_files,
    BPPARAM = BiocParallel::MulticoreParam(workers = n_cores, progressbar = TRUE),
    function(isc_path) {
      yaml_path <- paste0(isc_path, ".yaml")
      tryCatch(
        process_label_transfer_dataset(isc_path, yaml_path, out_lt_dir),
        error = function(e) {
          warning("Error processing ", basename(isc_path), ": ", conditionMessage(e))
          return(NULL)
        }
      )
    }
  )
  
  # Filter successful conversions
  results <- Filter(Negate(is.null), results)
  
  message_time("=".repeat(70))
  message_time(paste("✓ Successfully processed", length(results), "datasets"))
  message_time("=".repeat(70))
  message_time("")
  message_time("Label transfer splits saved to: data/processed/label_transfer/")
  message_time("")
  message_time("Next steps:")
  message_time("  1. Review splits: ls -lh data/processed/label_transfer/")
  message_time("  2. Run label transfer benchmark:")
  message_time("     cd label_transfer_task && Rscript -e 'targets::tar_make()'")
  message_time("")
  
  invisible(results)
}

if (!interactive()) {
  main()
}
