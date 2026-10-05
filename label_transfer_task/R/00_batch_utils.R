# label_transfer_task/R/00_batch_utils.R --- Batch label-transfer utilities
# One reference dataset vs. all other datasets (same annotation) combined as query.

source("R/00_utils.R")
source("R/00_between_utils.R")

lt_batch_query_label <- "ALL_OTHER_BATCHES"

lt_batch_data_processed_dir <- function() {
  ensure_dir(proj_path("data/processed/label_transfer_batch/pairs"))
}

lt_batch_results_root <- function() {
  ensure_dir(proj_path("label_transfer_task/results_batch"))
}

lt_batch_raw_results_dir <- function() {
  ensure_dir(file.path(lt_batch_results_root(), "raw_results"))
}

lt_batch_aggregated_dir <- function() {
  ensure_dir(file.path(lt_batch_results_root(), "aggregated"))
}

lt_batch_figures_dir <- function() {
  ensure_dir(file.path(lt_batch_results_root(), "figures"))
}

lt_batch_consistency_dir <- function() {
  ensure_dir(file.path(lt_batch_results_root(), "consistency"))
}

# One row per reference: pair_id, reference_dataset_id, query_dataset_id (label),
# query_dataset_ids (";"-separated members of the combined query), framework info.
list_batch_pairs <- function(params) {
  between_params <- params
  between_params$pairs_filter <- NULL
  bp <- list_between_dataset_pairs(between_params)

  empty <- tibble::tibble(
    pair_id = character(0),
    reference_dataset_id = character(0),
    query_dataset_id = character(0),
    query_dataset_ids = character(0),
    annotation_framework = character(0),
    annotation_reference = character(0),
    condition = character(0)
  )
  if (nrow(bp) == 0) return(empty)

  keys <- unique(bp[, c("reference_dataset_id", "annotation_framework", "annotation_reference")])
  out <- lapply(seq_len(nrow(keys)), function(i) {
    sub <- bp[bp$reference_dataset_id == keys$reference_dataset_id[i] &
                bp$annotation_framework == keys$annotation_framework[i], , drop = FALSE]
    tibble::tibble(
      pair_id = paste(keys$reference_dataset_id[i], "TO", lt_batch_query_label, sep = "__"),
      reference_dataset_id = keys$reference_dataset_id[i],
      query_dataset_id = lt_batch_query_label,
      query_dataset_ids = paste(sort(unique(sub$query_dataset_id)), collapse = ";"),
      annotation_framework = keys$annotation_framework[i],
      annotation_reference = keys$annotation_reference[i],
      condition = sub$condition[1]
    )
  })
  pairs <- dplyr::bind_rows(out)

  rf <- params$references_filter
  if (!is.null(rf) && length(rf) > 0) {
    pairs <- pairs[pairs$reference_dataset_id %in% as.character(rf), , drop = FALSE]
  }
  pairs
}

lt_batch_read_counts <- function(obj) {
  assay <- SeuratObject::DefaultAssay(obj)
  tryCatch(
    SeuratObject::GetAssayData(obj, assay = assay, layer = "counts"),
    error = function(e) SeuratObject::GetAssayData(obj, assay = assay, slot = "counts")
  )
}

lt_batch_load_dataset <- function(dataset_id) {
  path <- file.path(lt_isc_processed_dir(), paste0(dataset_id, ".rds"))
  if (!file.exists(path)) return(NULL)
  obj <- readRDS(path)
  if (!inherits(obj, "Seurat") || !"sample" %in% colnames(obj@meta.data)) return(NULL)
  md <- obj@meta.data
  ident <- lt_get_ident_for_dataset(dataset_id, colnames(md))
  list(
    counts = lt_batch_read_counts(obj),
    metadata = data.frame(
      sample = md$sample,
      cell_type = md[[ident]],
      dataset_id = dataset_id,
      stringsAsFactors = FALSE,
      row.names = rownames(md)
    )
  )
}

prepare_batch_pair <- function(pair_id,
                               reference_dataset_id,
                               query_dataset_ids,
                               query_dataset_id = lt_batch_query_label) {
  ref <- lt_batch_load_dataset(reference_dataset_id)
  if (is.null(ref)) {
    warning("Missing/invalid reference for ", pair_id)
    return(invisible(NULL))
  }

  q_ids <- strsplit(query_dataset_ids, ";", fixed = TRUE)[[1]]
  qs <- Filter(Negate(is.null), lapply(q_ids, lt_batch_load_dataset))
  if (length(qs) == 0) {
    warning("No valid query datasets for ", pair_id)
    return(invisible(NULL))
  }

  genes <- Reduce(intersect, lapply(qs, function(q) rownames(q$counts)))
  cell_names <- unlist(lapply(qs, function(q) colnames(q$counts)))
  prefix_cells <- anyDuplicated(cell_names) > 0 || any(cell_names %in% colnames(ref$counts))

  counts_list <- lapply(qs, function(q) {
    m <- q$counts[genes, , drop = FALSE]
    if (prefix_cells) colnames(m) <- paste(q$metadata$dataset_id[1], colnames(m), sep = "|")
    m
  })
  md_list <- lapply(qs, function(q) {
    md <- q$metadata
    md$source_dataset_id <- md$dataset_id
    md$dataset_id <- query_dataset_id
    if (prefix_cells) rownames(md) <- paste(md$source_dataset_id[1], rownames(md), sep = "|")
    md
  })

  query <- list(
    counts = do.call(cbind, counts_list),
    metadata = do.call(rbind, md_list)
  )
  ref$metadata$dataset_id <- reference_dataset_id

  outdir <- ensure_dir(file.path(lt_batch_data_processed_dir(), pair_id))
  saveRDS(ref, file.path(outdir, "reference.rds"))
  saveRDS(query, file.path(outdir, "query.rds"))

  invisible(file.path(outdir, "query.rds"))
}
