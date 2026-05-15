# label_transfer_task/R/00_between_utils.R --- Between-dataset label-transfer utilities

source("R/00_utils.R")

# ============================================================================
# BETWEEN-DATASET PATHS
# ============================================================================

lt_between_data_processed_dir <- function() {
  ensure_dir(proj_path("data/processed/label_transfer_between/pairs"))
}

lt_between_results_root <- function() {
  ensure_dir(proj_path("label_transfer_task/results_between"))
}

lt_between_raw_results_dir <- function() {
  ensure_dir(file.path(lt_between_results_root(), "raw_results"))
}

lt_between_aggregated_dir <- function() {
  ensure_dir(file.path(lt_between_results_root(), "aggregated"))
}

lt_between_figures_dir <- function() {
  ensure_dir(file.path(lt_between_results_root(), "figures"))
}

lt_between_consistency_dir <- function() {
  ensure_dir(file.path(lt_between_results_root(), "consistency"))
}

# ============================================================================
# BETWEEN-DATASET CONFIG
# ============================================================================

load_lt_between_params <- function() {
  load_pipeline_config("label_transfer_between_parameters.yaml")
}

get_lt_between_pairs_filter <- function(params) {
  pf <- params$pairs_filter
  if (is.null(pf) || length(pf) == 0) {
    return(NULL)
  }
  as.character(pf)
}

# ============================================================================
# SPECS MAPPING AND PAIRING
# ============================================================================

normalize_specs_token <- function(x) {
  x <- as.character(x)
  x <- gsub("[^a-zA-Z0-9._-]", ".", x)
  x <- gsub("\\.+", ".", x)
  x
}

prefix_from_annotation_reference <- function(annotation_reference) {
  ar <- as.character(annotation_reference)

  if (grepl("Stephenson", ar, fixed = TRUE)) return("Stephenson")
  if (grepl("Joanito", ar, fixed = TRUE)) return("JoaI")
  if (grepl("Gondal", ar, fixed = TRUE)) return("ICBAtlas")
  if (grepl("Sikkema", ar, fixed = TRUE)) return("LungAtlas")
  if (grepl("Andreatta", ar, fixed = TRUE)) return("BCC")

  NA_character_
}

specs_row_to_dataset_id <- function(row_df) {
  prefix <- prefix_from_annotation_reference(row_df[["Annotation reference"]])
  if (is.na(prefix)) {
    return(NA_character_)
  }

  batch <- normalize_specs_token(row_df[["Batch"]])
  condition <- normalize_specs_token(row_df[["Condition"]])

  if (prefix == "BCC") {
    return(paste(prefix, batch, "all", sep = "_"))
  }

  if (is.na(condition) || !nzchar(condition)) {
    return(paste(prefix, batch, sep = "_"))
  }

  paste(prefix, batch, condition, sep = "_")
}

list_between_dataset_pairs <- function(params = NULL) {
  if (is.null(params)) {
    params <- load_lt_between_params()
  }

  specs_path <- proj_path("data_processing/config/specs_datasets.csv")
  if (!file.exists(specs_path)) {
    stop("Missing specs file: ", specs_path)
  }

  specs <- read.csv(specs_path, stringsAsFactors = FALSE, check.names = FALSE)

  required_cols <- c("Annotation reference", "Condition", "Batch", "Label-Transfer Task")
  missing_cols <- setdiff(required_cols, colnames(specs))
  if (length(missing_cols) > 0) {
    stop("Missing columns in specs_datasets.csv: ", paste(missing_cols, collapse = ", "))
  }

  specs <- specs[specs[["Label-Transfer Task"]] == "yes", , drop = FALSE]
  specs$dataset_id <- vapply(seq_len(nrow(specs)), function(i) {
    specs_row_to_dataset_id(specs[i, , drop = FALSE])
  }, character(1))

  specs <- specs[!is.na(specs$dataset_id), , drop = FALSE]

  available_ids <- tools::file_path_sans_ext(
    list.files(lt_isc_processed_dir(), pattern = "\\.rds$", full.names = FALSE)
  )
  specs <- specs[specs$dataset_id %in% available_ids, , drop = FALSE]

  group_key <- paste(specs[["Annotation reference"]], specs[["Condition"]], sep = "__")
  groups <- split(specs, group_key)

  all_pairs <- lapply(groups, function(g) {
    if (nrow(g) < 2) return(NULL)

    g <- g[!duplicated(g$dataset_id), , drop = FALSE]
    if (nrow(g) < 2) return(NULL)

    cmb <- utils::combn(seq_len(nrow(g)), 2)

    # Build ordered pairs so each dataset can be reference and query.
    fwd <- data.frame(
      reference_dataset_id = g$dataset_id[cmb[1, ]],
      query_dataset_id = g$dataset_id[cmb[2, ]],
      reference_batch = g[["Batch"]][cmb[1, ]],
      query_batch = g[["Batch"]][cmb[2, ]],
      annotation_reference = g[["Annotation reference"]][cmb[1, ]],
      condition = g[["Condition"]][cmb[1, ]],
      stringsAsFactors = FALSE
    )
    rev <- data.frame(
      reference_dataset_id = g$dataset_id[cmb[2, ]],
      query_dataset_id = g$dataset_id[cmb[1, ]],
      reference_batch = g[["Batch"]][cmb[2, ]],
      query_batch = g[["Batch"]][cmb[1, ]],
      annotation_reference = g[["Annotation reference"]][cmb[1, ]],
      condition = g[["Condition"]][cmb[1, ]],
      stringsAsFactors = FALSE
    )

    out <- rbind(fwd, rev)
    out[out$reference_batch != out$query_batch, , drop = FALSE]
  })

  pairs <- do.call(rbind, all_pairs)
  if (is.null(pairs) || nrow(pairs) == 0) {
    return(tibble::tibble(
      pair_id = character(0),
      reference_dataset_id = character(0),
      query_dataset_id = character(0),
      annotation_reference = character(0),
      condition = character(0)
    ))
  }

  pairs$pair_id <- paste(pairs$reference_dataset_id, "TO", pairs$query_dataset_id, sep = "__")

  pair_filter <- get_lt_between_pairs_filter(params)
  if (!is.null(pair_filter)) {
    pairs <- pairs[pairs$pair_id %in% pair_filter, , drop = FALSE]
  }

  tibble::as_tibble(unique(pairs[, c(
    "pair_id",
    "reference_dataset_id",
    "query_dataset_id",
    "annotation_reference",
    "condition"
  )]))
}

# ============================================================================
# PREPARE BETWEEN-DATASET REFERENCE/QUERY FILES
# ============================================================================

prepare_between_dataset_pair <- function(pair_id,
                                         reference_dataset_id,
                                         query_dataset_id,
                                         params = NULL) {
  if (is.null(params)) {
    params <- load_lt_between_params()
  }

  ref_path <- file.path(lt_isc_processed_dir(), paste0(reference_dataset_id, ".rds"))
  qry_path <- file.path(lt_isc_processed_dir(), paste0(query_dataset_id, ".rds"))

  if (!file.exists(ref_path) || !file.exists(qry_path)) {
    warning("Missing pair inputs for ", pair_id)
    return(invisible(NULL))
  }

  ref_obj <- readRDS(ref_path)
  qry_obj <- readRDS(qry_path)

  if (!inherits(ref_obj, "Seurat") || !inherits(qry_obj, "Seurat")) {
    warning("Pair contains non-Seurat object: ", pair_id)
    return(invisible(NULL))
  }

  ref_md <- ref_obj@meta.data
  qry_md <- qry_obj@meta.data

  if (!"sample" %in% colnames(ref_md) || !"sample" %in% colnames(qry_md)) {
    warning("Missing 'sample' metadata in pair: ", pair_id)
    return(invisible(NULL))
  }

  ref_ident <- lt_get_ident_for_dataset(reference_dataset_id, colnames(ref_md))
  qry_ident <- lt_get_ident_for_dataset(query_dataset_id, colnames(qry_md))

  ref_counts <- tryCatch(
    SeuratObject::GetAssayData(
      ref_obj,
      assay = SeuratObject::DefaultAssay(ref_obj),
      layer = "counts"
    ),
    error = function(e) {
      SeuratObject::GetAssayData(
        ref_obj,
        assay = SeuratObject::DefaultAssay(ref_obj),
        slot = "counts"
      )
    }
  )

  qry_counts <- tryCatch(
    SeuratObject::GetAssayData(
      qry_obj,
      assay = SeuratObject::DefaultAssay(qry_obj),
      layer = "counts"
    ),
    error = function(e) {
      SeuratObject::GetAssayData(
        qry_obj,
        assay = SeuratObject::DefaultAssay(qry_obj),
        slot = "counts"
      )
    }
  )

  reference <- list(
    counts = ref_counts,
    metadata = data.frame(
      sample = ref_md$sample,
      cell_type = ref_md[[ref_ident]],
      dataset_id = reference_dataset_id,
      stringsAsFactors = FALSE,
      row.names = rownames(ref_md)
    )
  )

  query <- list(
    counts = qry_counts,
    metadata = data.frame(
      sample = qry_md$sample,
      cell_type = qry_md[[qry_ident]],
      dataset_id = query_dataset_id,
      stringsAsFactors = FALSE,
      row.names = rownames(qry_md)
    )
  )

  outdir <- ensure_dir(file.path(lt_between_data_processed_dir(), pair_id))
  saveRDS(reference, file.path(outdir, "reference.rds"))
  saveRDS(query, file.path(outdir, "query.rds"))

  invisible(file.path(outdir, "query.rds"))
}
