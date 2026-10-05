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
# The query contains every other eligible dataset with the same annotation,
# including datasets from the reference's own batch (e.g. Bassez Pre for Bassez Post).
list_batch_pairs <- function(params) {
  specs_path <- proj_path("data_processing/config/specs_datasets.csv")
  if (!file.exists(specs_path)) stop("Missing specs file: ", specs_path)
  specs <- read.csv(specs_path, stringsAsFactors = FALSE, check.names = FALSE)

  specs <- specs[specs[["Label-Transfer Task"]] == "yes", , drop = FALSE]
  specs$dataset_id <- vapply(seq_len(nrow(specs)), function(i) {
    specs_row_to_dataset_id(specs[i, , drop = FALSE])
  }, character(1))
  specs <- specs[!is.na(specs$dataset_id) & !duplicated(specs$dataset_id), , drop = FALSE]

  available_ids <- tools::file_path_sans_ext(
    list.files(lt_isc_processed_dir(), pattern = "\\.rds$", full.names = FALSE)
  )
  specs <- specs[specs$dataset_id %in% available_ids, , drop = FALSE]

  groups <- split(specs, paste(specs[["# Annotation frameworks"]], specs[["Annotation reference"]], sep = "||"))
  out <- lapply(groups, function(g) {
    if (nrow(g) < 2) return(NULL)
    dplyr::bind_rows(lapply(seq_len(nrow(g)), function(i) {
      tibble::tibble(
        pair_id = paste(g$dataset_id[i], "TO", lt_batch_query_label, sep = "__"),
        reference_dataset_id = g$dataset_id[i],
        query_dataset_id = lt_batch_query_label,
        query_dataset_ids = paste(sort(g$dataset_id[-i]), collapse = ";"),
        annotation_framework = g[["# Annotation frameworks"]][i],
        annotation_reference = g[["Annotation reference"]][i],
        condition = g[["Condition"]][i]
      )
    }))
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
