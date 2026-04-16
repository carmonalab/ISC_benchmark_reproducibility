# label_transfer_task/R/00_utils.R --- Label-transfer-specific utility functions

source("../R/shared_helpers.R")

# ============================================================================
# LABEL-TRANSFER-SPECIFIC PATHS
# ============================================================================

lt_data_processed_dir <- function(replicate = NULL) {
  base <- proj_path("data/processed/label_transfer")
  if (is.null(replicate)) {
    return(base)
  }
  ensure_dir(file.path(base, paste0("rep", replicate)))
}

lt_isc_processed_dir <- function() {
  proj_path("data/processed")
}

lt_results_root <- function() {
  ensure_dir(proj_path("label_transfer_task/results"))
  proj_path("label_transfer_task/results")
}

lt_raw_results_dir <- function() {
  ensure_dir(file.path(lt_results_root(), "raw_results"))
}

lt_aggregated_dir <- function() {
  ensure_dir(file.path(lt_results_root(), "aggregated"))
}

lt_figures_dir <- function() {
  ensure_dir(file.path(lt_results_root(), "figures"))
}

lt_consistency_dir <- function() {
  ensure_dir(file.path(lt_results_root(), "consistency"))
}

# ============================================================================
# LOAD LABEL-TRANSFER PARAMETERS
# ============================================================================

load_lt_params <- function() {
  # Load all label transfer parameters from pipeline-specific config
  # (All required settings are self-contained in label_transfer_parameters.yaml)
  load_pipeline_config("label_transfer_parameters.yaml")
}

get_lt_classifiers <- function(params) {
  # Return list of classifiers to evaluate
  if (is.null(params$classifiers$methods) || length(params$classifiers$methods) == 0) {
    stop("Missing required config: classifiers.methods")
  }
  params$classifiers$methods
}

get_lt_n_replicates <- function(params) {
  n <- params$n_replicates
  if (is.null(n) || !is.numeric(n) || length(n) != 1 || is.na(n) || n < 1) {
    stop("Missing/invalid required config: n_replicates")
  }
  as.integer(n)
}

get_lt_seed <- function(params) {
  s <- params$seed
  if (is.null(s) || !is.numeric(s) || length(s) != 1 || is.na(s)) {
    stop("Missing/invalid required config: seed")
  }
  as.integer(s)
}

get_lt_n_cores <- function(params) {
  n <- params$n_cores
  if (is.null(n) || !is.numeric(n) || length(n) != 1 || is.na(n) || n < 1) {
    stop("Missing/invalid required config: n_cores")
  }
  as.integer(n)
}

get_lt_data_prep <- function(params) {
  dp <- params$data_preparation
  if (is.null(dp$min_samples) || is.null(dp$prop_query)) {
    stop("Missing required config: data_preparation.min_samples and/or data_preparation.prop_query")
  }
  if (!is.numeric(dp$min_samples) || dp$min_samples < 1) {
    stop("Invalid config: data_preparation.min_samples")
  }
  if (!is.numeric(dp$prop_query) || dp$prop_query <= 0 || dp$prop_query >= 1) {
    stop("Invalid config: data_preparation.prop_query (must be between 0 and 1)")
  }
  dp
}

# ============================================================================
# DATASET IDENTIFICATION COLUMN MAPPING
# ============================================================================

#' Get cell-type identification column name by dataset prefix
#'
#' Maps dataset identifiers to their corresponding cell-type annotation columns.
#' This ensures consistent identification across different datasets.
#'
#' @param prefix Character: Dataset prefix (e.g., "JoaI", "Mitchel", "BCC")
#'
#' @return Character: Column name to use for cell-type identification
#'
#' @details
#' Maps dataset prefixes to their metadata columns:
#' - JoaI (Joanito et al., 2022 CRC): "cell.type"
#' - Mitchel (Mitchel et al., 2023): "OriginalAnnotationLevel1"
#' - BCC (Yost et al., BCC atlas): "annotation"
#' - LungAtlas (Sikkema et al., HCA Lung): "cell_type"
#' - ICBAtlas (Gondal et al., ICB Atlas): "cell_type"
#'
get_idents_by_prefix <- function(prefix) {
  switch(prefix,
    "JoaI"        = "cell.type",
    "StephensonE" = "OriginalAnnotationLevel1",
    "BCC"         = "annotation",
    "LungAtlas"   = "cell_type",
    "ICBAtlas"    = "cell_type",
    "celltype"    # default fallback
  )
}

# ============================================================================
# VALIDATION & DATA PREPARATION
# ============================================================================

validate_label_transfer_participation <- function(dataset_id) {
  # Check if dataset participates in label_transfer task
  meta_path <- file.path(
    proj_path("data/processed"),
    paste0(dataset_id, ".rds.yaml")
  )
  
  if (!file.exists(meta_path)) {
    return(TRUE)  # Default: include if metadata missing
  }
  
  meta <- load_config(meta_path)
  isTRUE(meta$label_transfer)
}

validate_label_transfer_data <- function(dataset_id) {
  # Check if query + reference files exist
  data_dir <- lt_data_processed_dir()
  query_path <- file.path(data_dir, dataset_id, "query.rds")
  ref_path <- file.path(data_dir, dataset_id, "reference.rds")
  
  all(file.exists(query_path, ref_path))
}

validate_label_transfer_data_replicate <- function(dataset_id, replicate) {
  data_dir <- lt_data_processed_dir(replicate)
  query_path <- file.path(data_dir, dataset_id, "query.rds")
  ref_path <- file.path(data_dir, dataset_id, "reference.rds")
  all(file.exists(query_path, ref_path))
}

list_label_transfer_datasets_from_isc <- function(params = NULL) {
  if (is.null(params)) {
    params <- load_lt_params()
  }

  dp <- get_lt_data_prep(params)
  min_samples <- as.integer(dp$min_samples)

  isc_dir <- lt_isc_processed_dir()
  if (!dir.exists(isc_dir)) {
    warning("ISC processed directory not found: ", isc_dir)
    return(character(0))
  }

  isc_files <- list.files(isc_dir, pattern = "\\.rds$", full.names = FALSE, ignore.case = TRUE)
  dataset_ids <- tools::file_path_sans_ext(isc_files)

  dataset_ids <- dataset_ids[vapply(dataset_ids, validate_label_transfer_participation, logical(1))]

  # Filter by sample threshold (so downstream split-prep always produces files)
  keep <- vapply(dataset_ids, function(id) {
    path <- file.path(isc_dir, paste0(id, ".rds"))
    obj <- tryCatch(readRDS(path), error = function(e) NULL)
    if (is.null(obj) || !inherits(obj, "Seurat")) {
      return(FALSE)
    }
    md <- obj@meta.data
    if (!"sample" %in% colnames(md)) {
      return(FALSE)
    }
    length(unique(md$sample)) >= min_samples
  }, logical(1))

  dataset_ids <- dataset_ids[keep]

  # Optional filter from config
  filter <- params$datasets_filter
  if (!is.null(filter) && length(filter) > 0) {
    dataset_ids <- dataset_ids[dataset_ids %in% filter]
  }

  dataset_ids
}

lt_get_ident_for_dataset <- function(dataset_key, metadata_cols) {
  prefix <- strsplit(dataset_key, "_")[[1]][1]
  ident_col <- get_idents_by_prefix(prefix)

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

prepare_label_transfer_split <- function(dataset_id, replicate, params) {
  dp <- get_lt_data_prep(params)
  min_samples <- as.integer(dp$min_samples)
  prop_query <- as.numeric(dp$prop_query)

  seed_base <- get_lt_seed(params)
  split_seed <- seed_base + as.integer(replicate) - 1L

  isc_path <- file.path(lt_isc_processed_dir(), paste0(dataset_id, ".rds"))
  if (!file.exists(isc_path)) {
    stop("ISC object not found: ", isc_path)
  }

  yaml_path <- paste0(isc_path, ".yaml")
  if (file.exists(yaml_path)) {
    meta <- yaml::read_yaml(yaml_path)
    if (!isTRUE(meta$label_transfer)) {
      return(invisible(NULL))
    }
  }

  obj <- readRDS(isc_path)
  if (!inherits(obj, "Seurat")) {
    warning("Not a Seurat object: ", isc_path)
    return(invisible(NULL))
  }

  if (!requireNamespace("SeuratObject", quietly = TRUE)) {
    stop("Package 'SeuratObject' is required to extract counts")
  }

  counts <- SeuratObject::GetAssayData(
    obj,
    assay = SeuratObject::DefaultAssay(obj),
    slot = "counts"
  )

  md <- obj@meta.data
  if (!"sample" %in% colnames(md)) {
    stop("Missing required metadata column 'sample' in ", dataset_id)
  }

  n_samples <- length(unique(md$sample))
  if (n_samples < min_samples) {
    return(invisible(NULL))
  }

  ident_col <- lt_get_ident_for_dataset(dataset_id, colnames(md))

  md_standardized <- data.frame(
    sample = md$sample,
    cell_type = md[[ident_col]],
    dataset_id = dataset_id,
    row.names = rownames(md),
    stringsAsFactors = FALSE
  )

  all_samples <- unique(md_standardized$sample)
  n_query <- floor(length(all_samples) * prop_query)

  set.seed(split_seed)
  query_samples <- sample(all_samples, size = n_query)

  cells_query <- rownames(md_standardized)[md_standardized$sample %in% query_samples]
  cells_ref <- rownames(md_standardized)[!md_standardized$sample %in% query_samples]

  out_base <- lt_data_processed_dir(replicate)
  outdir <- ensure_dir(file.path(out_base, dataset_id))

  query <- list(
    counts = counts[, cells_query, drop = FALSE],
    metadata = md_standardized[cells_query, , drop = FALSE]
  )
  reference <- list(
    counts = counts[, cells_ref, drop = FALSE],
    metadata = md_standardized[cells_ref, , drop = FALSE]
  )

  saveRDS(query, file.path(outdir, "query.rds"))
  saveRDS(reference, file.path(outdir, "reference.rds"))

  invisible(file.path(outdir, "query.rds"))
}

# ============================================================================
# HELPER: LOAD CLASSIFIER DEFINITIONS
# ============================================================================

load_classifier_code <- function() {
  # Source the classifiers.R file that defines classifier functions
  classifier_file <- proj_path("label_transfer_task/classifiers/classifiers.R")
  
  if (!file.exists(classifier_file)) {
    stop("Classifier definitions not found: ", classifier_file)
  }
  
  source(classifier_file, local = TRUE)
  invisible(NULL)
}
