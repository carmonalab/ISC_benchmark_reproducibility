suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(Seurat)
  library(SeuratObject)
  library(BiocParallel)
})

# Source shared utilities (defines proj_path and normalize_metadata_name)
source("../R/shared_helpers.R")
source(proj_path("R/cli_utils.R"))


#' Load dataset from raw file(s)

load_dataset <- function(ds, input_dir = "../data/raw") {
  tryCatch({
    # Single file (default)
    if (!is.null(ds$raw_file)) {
      file <- file.path(input_dir, ds$raw_file)
      return(readRDS(file))
    }
    
    return(NULL)
  }, error = function(e) {
    warning(sprintf("Failed to load dataset '%s': %s", ds$id, e$message))
    return(NULL)
  })
}

#' Process one dataset: downsample → split by batch+condition → save
#'
#' @param obj Seurat object
#' @param prefix Output file prefix (e.g., "JoaI", "Mitchel")
#' @param ident_col Cell type annotation column
#' @param sample_col Sample identifier column
#' @param batch_col Batch/site/dataset identifier column
#' @param condition_col Condition column (e.g., status, tissue, pre_post)
#' @param exclude_celltypes Optional: cell types to exclude
#' @param config Configuration list with reproducibility and filtering parameters
#'
#' @return NULL (saves files directly)
process_dataset <- function(obj,
                           prefix,
                           ident_col,
                           sample_col,
                           batch_col = NULL,
                           condition_col = NULL,
                           exclude_celltypes = NULL,
                           config) {
  
  message(sprintf("[%s] Processing %s - %d cells", 
                  format(Sys.time(), "%H:%M:%S"), prefix, ncol(obj)))
  
  # ===== Step 1: Exclude low consistency cell types =====
  if (!is.null(exclude_celltypes) && length(exclude_celltypes) > 0) {
    if (ident_col %in% colnames(obj@meta.data)) {
      n_before <- ncol(obj)
      keep <- !obj[[ident_col, drop = TRUE]] %in% exclude_celltypes
      obj <- obj[, keep]
      n_removed <- n_before - ncol(obj)
      if (n_removed > 0) {
        message(sprintf("  Removed %d cells of excluded types: %s", 
                       n_removed, paste(exclude_celltypes, collapse = ", ")))
      }
    }
  }
  
  # ===== Step 2: Standardize sample column =====
  obj$sample <- obj[[sample_col, drop = TRUE]]
  
  # ===== Step 3: Downsample by cell type + sample =====
  message(sprintf("  Downsampling (max %d per cell type per sample)", 
                  config$downsampling$max_cells_per_sample_celltype))
  obj$Sample_annot <- paste(obj$sample, obj[[ident_col, drop = TRUE]])
  Idents(obj) <- "Sample_annot"
  
  set.seed(config$seed)
  obj <- subset(obj, downsample = config$downsampling$max_cells_per_sample_celltype)
  
  message(sprintf("  After downsampling: %d cells", ncol(obj)))
  
  # ===== Step 4: Create split variable (batch_condition) =====

has_batch <- !is.null(batch_col) && batch_col %in% colnames(obj@meta.data)
has_condition <- !is.null(condition_col) && condition_col %in% colnames(obj@meta.data)

if (has_batch && has_condition) {
  # Both batch and condition
  batch_vals <- normalize_metadata_name(obj[[batch_col, drop = TRUE]])
  condition_vals <- normalize_metadata_name(obj[[condition_col, drop = TRUE]])
  obj$split <- paste(batch_vals, condition_vals, sep = "_")
  
} else if (has_batch) {
  # Batch only
  obj$split <- normalize_metadata_name(obj[[batch_col, drop = TRUE]])
  
} else if (has_condition) {
  # Condition only
  obj$split <- normalize_metadata_name(obj[[condition_col, drop = TRUE]])
  
} else {
  # Neither provided → fallback
  obj$split <- "all"
}
  
  # ===== Step 5: Split by batch+condition =====
  message("  Splitting by batch+condition")
  splits <- SplitObject(obj, split.by = "split")
  message(sprintf("  Created %d splits", length(splits)))
  
  # ===== Step 6: Apply optimal_dataset to each split =====
  message("  Applying optimal_dataset filtering")
  
  w <- min(config$n_cores, length(splits))
  # if resulting splits are too large we prioritize keeping samples
  # with enough cell types and total cells
  splits <- bplapply(
    splits,
    BPPARAM = MulticoreParam(workers = w, progressbar = FALSE),
    function(x) {
      optimal_dataset(
        x,
        ident = ident_col,
        sample = "sample",
        min_samples = config$sample_filtering$min_samples,
        max_samples = config$sample_filtering$max_samples,
        min_cells_sample = config$sample_filtering$min_cells_per_sample
      )
    }
  )
  
  # Filter NULL results (failed optimal_dataset checks)
  splits <- Filter(Negate(is.null), splits)
  
  if (length(splits) == 0) {
    warning("No valid splits after optimal_dataset filtering for ", prefix)
    return(invisible(NULL))
  }
  
  message(sprintf("  Valid splits after filtering: %d", length(splits)))
  
  # ===== Step 7: Save each split =====
  message("  Saving split files")
  
  dir.create(config$out_dir, showWarnings = FALSE, recursive = TRUE)
  
  w <- min(config$n_cores, length(splits))
  bplapply(
    names(splits),
    BPPARAM = MulticoreParam(workers = w, progressbar = FALSE),
    function(name) {
      out_file <- file.path(config$out_dir, paste0(prefix, "_", name, ".rds"))
      saveRDS(splits[[name]], out_file)
    }
  )
  
  message(sprintf("✓ Saved %d processed datasets for %s", length(splits), prefix))
  invisible(NULL)
}

#' Filter Seurat object to enforce minimum sample/cell requirements
#'
#' Applies quality control filters to ensure sufficient samples and cell representation.
#' This function performs two main checks:
#' 1. **Sample count validation**: Ensures the object has at least `min_samples` samples.
#'    If exceeds `max_samples`, selects top samples by cell type diversity and total cells.
#' 2. **Cell representation**: Preferentially retains samples with:
#'    - Multiple cell types (nident > 1)
#'    - Sufficient total cells (>= min_cells_sample)
#'
#' This is commonly used to filter batch/condition-specific splits after downsampling,
#' ensuring each split has balanced representation across samples and cell types.
#'
#' @param obj Seurat object to filter
#' @param ident Character: Column name in metadata containing cell type/identity annotations
#'   (default: "celltype")
#' @param sample Character: Column name in metadata containing sample identifiers
#'   (default: "sample")
#' @param min_samples Numeric: Minimum required samples. Returns NULL if fewer samples
#'   remain after filtering (default: 9)
#' @param max_samples Numeric: Maximum samples to retain. If exceeded, selects top samples
#'   by cell type diversity and cell count (default: 15)
#' @param min_cells_sample Numeric: Minimum cells required per sample to be eligible for
#'   selection when subsampling (default: 200)
#' @param verbose Logical: If TRUE, print filtering messages (default: TRUE)
#'
#' @return
#'   - Seurat object (filtered) if sample count passes min/max checks
#'   - NULL if final sample count < min_samples (indicating split is not suitable)
#'
#' @details
#' **Filtering Algorithm:**
#' 1. Count unique samples and cell types per sample
#' 2. If samples > max_samples: select top `max_samples` by:
#'    - Priority 1: Cell type diversity (samples with multiple cell types)
#'    - Priority 2: Total cell count (samples with more cells)
#'    - Filter: Exclude samples with < min_cells_sample cells
#' 3. Subset object to selected samples
#' 4. Check final sample count:
#'    - If < min_samples: return NULL (split discarded)
#'    - If >= min_samples: return filtered object
#'
#' @examples
#' \dontrun{
#' # Filter a dataset split to ensure balanced representation
#' filtered_obj <- optimal_dataset(
#'   obj = split_obj,
#'   ident = "cell_type",
#'   sample = "sample",
#'   min_samples = 9,
#'   max_samples = 15,
#'   min_cells_sample = 200
#' )
#' }
optimal_dataset <- function(obj,
                            ident = "celltype",
                            sample = "sample",
                            min_samples = 9,
                            max_samples = 15,
                            min_cells_sample = 200,
                            verbose = TRUE) {
  md <- obj@meta.data %>%
    mutate(ident = .data[[ident]],
           sample = .data[[sample]]) %>%
    select(sample, ident)

  nsamples <- length(unique(md[["sample"]]))

  nidents <- md %>%
    group_by(sample, ident) %>%
    summarize(ncells_ident = n(), .groups = "drop") %>%
    group_by(sample) %>%
    mutate(nident = n(),
           ncells_sample = sum(ncells_ident))

  if (nsamples > max_samples) {
    top_rep <- nidents %>%
      distinct(sample, nident, ncells_sample) %>%
      filter(nident > 1 & ncells_sample >= min_cells_sample) %>%
      arrange(desc(nident), desc(ncells_sample)) %>%
      head(max_samples) %>%
      pull(sample)

    obj <- obj[, obj@meta.data[[sample]] %in% top_rep]
    if (verbose) {
      message("- Subsetting dataset from ", nsamples, " to ", length(top_rep), " samples")
    }
    nsamples <- length(unique(obj@meta.data[[sample]]))
  }

  if (nsamples < min_samples) {
    if (verbose) message("- Only ", nsamples, " samples; returning NULL")
    return(NULL)
  }

  if (verbose) message("- Returning ", nsamples, " samples")
  obj
}


#' Apply dataset-specific preprocessing
#' @param obj Seurat object
#' @param ds Dataset configuration entry
#' @return Processed Seurat object
preprocess_dataset <- function(obj, ds) {
  
  # JoaI: Filter cancer subtypes and standardize condition
  if (ds$id == "JoaI_2022") {
    obj <- obj[, !obj$iCMS %in% c("iCMS2", "iCMS3")]
    obj$sample.origin <- ifelse(
      grepl("tumor", obj$sample.origin, ignore.case = TRUE),
      "Tumor", 
      "Normal"
    )
  }
  
  # LungAtlas: Clean tissue and dataset columns
  if (ds$id == "LungAtlas_2023") {
    obj$tissue <- normalize_metadata_name(obj$tissue)
    obj$dataset <- normalize_metadata_name(obj$dataset)
  }
  
  # ICBAtlas: Clean condition and dataset columns
  if (ds$id == "ICBAtlas_2024") {
    obj$pre_post <- normalize_metadata_name(obj$pre_post)
    
    # Clean the batch column used in core_datasets.yaml
    if ("Study_name" %in% colnames(obj@meta.data)) {
      obj$Study_name <- normalize_metadata_name(obj$Study_name)
    }
    # Backwards-compat: some objects may already have a `dataset` column
    if ("dataset" %in% colnames(obj@meta.data)) {
      obj$dataset <- normalize_metadata_name(obj$dataset)
    }
  }
  
  obj
}
