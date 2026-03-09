#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(yaml)
  library(dplyr)
  library(stringr)
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
})

source("utils/project_paths.R")
source("utils/cli_utils.R")
source("utils/data_processing_helpers.R")

stop_if_not_root <- function() {
  if (!file.exists("config/datasets.yaml") || !file.exists("data_processing/specs_datasets.csv")) {
    stop("Run this script from the repository root.")
  }
}

normalize_text <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  # Normalize various hidden/zero-width spaces often found in copy/pasted references.
  x <- stringr::str_replace_all(x, "[\u00A0\u200B\u200C\u200D\uFEFF]", " ")
  x <- stringr::str_replace_all(x, "[[:space:]]+", " ")
  x <- stringr::str_trim(x)
  x
}

normalize_match_token <- function(x) {
  x <- normalize_text(x)
  x <- tolower(x)
  x <- stringr::str_replace_all(x, "[^a-z0-9]+", "_")
  x <- stringr::str_replace_all(x, "^_+|_+$", "")
  x
}

match_meta_value <- function(meta_values, spec_value) {
  spec_value <- normalize_text(spec_value)
  if (!nzchar(spec_value)) return(rep(TRUE, length(meta_values)))

  spec_n <- normalize_match_token(spec_value)
  meta_n <- normalize_match_token(meta_values)
  if (!nzchar(spec_n)) return(rep(TRUE, length(meta_values)))

  keep <- meta_n == spec_n
  if (any(keep, na.rm = TRUE)) return(keep)

  # Fallback: prefix match either direction (handles "Banovich.Kropski" vs "Banovich_Kropski_2020").
  keep <- startsWith(meta_n, spec_n) | startsWith(spec_n, meta_n)
  keep
}

match_meta_any <- function(meta_values, spec_values) {
  spec_values <- as.character(spec_values)
  spec_values[is.na(spec_values)] <- ""
  spec_values <- normalize_text(spec_values)
  spec_values <- spec_values[nzchar(spec_values)]
  if (length(spec_values) == 0) return(rep(TRUE, length(meta_values)))

  keep <- rep(FALSE, length(meta_values))
  for (v in unique(spec_values)) {
    keep <- keep | match_meta_value(meta_values, v)
  }
  keep
}

sanitize_token <- function(x) {
  x <- normalize_text(x)
  x <- stringr::str_replace_all(x, "[[:space:]]+", ".")
  x <- stringr::str_replace_all(x, "[^A-Za-z0-9._-]", "")
  x
}

read_specs <- function(path) {
  read.csv(path, check.names = FALSE, stringsAsFactors = FALSE) %>%
    mutate(
      `ISC-benchmarking` = tolower(`ISC-benchmarking`),
      `Batch comparison` = tolower(`Batch comparison`),
      `Perturbation comparison` = tolower(`Perturbation comparison`),
      `Label-Transfer Task` = tolower(`Label-Transfer Task`)
    )
}

read_dataset_registry <- function(path) {
  cfg <- yaml::read_yaml(path)
  datasets <- cfg$datasets
  if (is.null(datasets) || length(datasets) == 0) stop("No datasets found in config/datasets.yaml")
  datasets
}

read_zenodo_doi <- function(path) {
  cfg <- yaml::read_yaml(path)
  doi <- cfg$zenodo$doi
  doi <- normalize_text(doi)
  if (!nzchar(doi) || identical(doi, "TODO")) return(NA_character_)
  doi
}

registry_by_reference <- function(datasets) {
  refs <- vapply(datasets, function(d) normalize_text(d$reference), character(1))
  if (anyDuplicated(refs)) {
    dup <- unique(refs[duplicated(refs)])
    stop(
      "Duplicate normalized references in config/datasets.yaml: ",
      paste(dup, collapse = "; "),
      "\nPlease disambiguate references (even if only whitespace differs)."
    )
  }
  stats::setNames(datasets, refs)
}

download_if_needed <- function(dataset_cfg, zenodo_raw_dir, download = FALSE) {
  if (isFALSE(download)) return(invisible(TRUE))

  files <- dataset_cfg$files
  if (is.null(files) || length(files) == 0) return(invisible(TRUE))

  for (fi in files) {
    if (is.null(fi$url) || identical(fi$url, "TODO")) {
      message("[download] Skipping (no url): ", fi$filename)
      next
    }
    dest <- file.path(zenodo_raw_dir, fi$filename)
    ensure_dir(dirname(dest))
    if (file.exists(dest)) {
      message("[download] Already present: ", dest)
      next
    }
    message("[download] Fetching: ", fi$filename)
    utils::download.file(fi$url, destfile = dest, mode = "wb", quiet = FALSE)
  }

  invisible(TRUE)
}

process_one_spec_group <- function(group_row,
                                  dataset_cfg,
                                  processing_cfg,
                                  label_transfer_cfg,
                                  zenodo_raw_dir,
                                  out_isc_dir,
                                  out_lt_dir,
                                  zenodo_doi = NA_character_) {
  used_any <- isTRUE(group_row$isc) || isTRUE(group_row$batch_comparison) ||
    isTRUE(group_row$perturbation_comparison) || isTRUE(group_row$label_transfer)
  if (!used_any) {
    message("Skipping (not used by any task): ", dataset_cfg$id, " / ", group_row$Batch, " / ", group_row$Condition)
    return(invisible(NULL))
  }

  raw_path <- file.path(zenodo_raw_dir, dataset_cfg$raw_filename)
  if (!file.exists(raw_path)) {
    doi_msg <- ""
    if (!is.na(zenodo_doi) && nzchar(zenodo_doi)) {
      doi_msg <- paste0("\nZenodo DOI: ", zenodo_doi)
    }
    stop(
      "Missing raw dataset file: ", raw_path,
      "\nSee data_processing/README.md and config/datasets.yaml",
      doi_msg
    )
  }

  message_time(paste("Loading raw dataset:", dataset_cfg$id))
  obj <- readRDS(raw_path)
  if (!inherits(obj, "Seurat")) {
    stop("Expected a Seurat object in ", raw_path)
  }

  md <- obj@meta.data
  for (col in c(dataset_cfg$batch_col, dataset_cfg$condition_col, dataset_cfg$sample_col)) {
    if (!is.null(col) && !col %in% colnames(md)) {
      stop("Column '", col, "' not found in metadata for dataset ", dataset_cfg$id)
    }
  }

  batch_vals <- group_row$Batch
  condition_vals <- group_row$Condition
  batch_vals <- as.character(batch_vals)
  condition_vals <- as.character(condition_vals)
  batch_vals[is.na(batch_vals)] <- ""
  condition_vals[is.na(condition_vals)] <- ""

  # Subset for batch/condition if specified in specs.
  if (!is.null(dataset_cfg$batch_col)) {
    keep <- match_meta_any(md[[dataset_cfg$batch_col]], batch_vals)
    obj <- obj[, keep]
  }
  if (!is.null(dataset_cfg$condition_col)) {
    md2 <- obj@meta.data
    keep <- match_meta_any(md2[[dataset_cfg$condition_col]], condition_vals)
    obj <- obj[, keep]
  }

  if (ncol(obj) == 0) {
    warning(
      "No cells after subsetting for ",
      dataset_cfg$id, " / ",
      paste(batch_vals, collapse = ";"), " / ",
      paste(condition_vals, collapse = ";")
    )
    return(invisible(NULL))
  }

  # Standardize commonly used metadata columns.
  obj$sample <- obj@meta.data[[dataset_cfg$sample_col]]
  obj$batch <- if (!is.null(dataset_cfg$batch_col)) obj@meta.data[[dataset_cfg$batch_col]] else NA
  obj$condition <- if (!is.null(dataset_cfg$condition_col)) obj@meta.data[[dataset_cfg$condition_col]] else NA
  obj$dataset_reference <- dataset_cfg$reference
  obj$dataset_id <- dataset_cfg$id
  # Used by batch-effect scripts in the original pipeline.
  obj$split <- obj$batch

  # Determine which annotation column to use for downsampling + label transfer.
  idents <- group_row$idents
  idents <- idents[idents != ""]
  ident_pref <- tail(idents, 1)
  ident_available <- idents[idents %in% colnames(obj@meta.data)]
  if (length(ident_available) == 0) {
    stop("None of the requested idents are present in metadata for ", dataset_cfg$id,
         ". Requested: ", paste(idents, collapse = ", "))
  }
  if (!ident_pref %in% ident_available) {
    ident_pref <- ident_available[[1]]
  }

  message_time(paste("Downsampling using ident:", ident_pref))
  obj$Sample_annot <- paste(obj$sample, obj@meta.data[[ident_pref]])
  Idents(obj) <- "Sample_annot"
  set.seed(processing_cfg$seed)
  obj2 <- subset(obj, downsample = processing_cfg$max_cells_per_sample_celltype)

  # Keep sample column name consistent.
  obj2$sample <- obj2@meta.data[["sample"]]

  if (isTRUE(processing_cfg$enforce_optimal_dataset)) {
    message_time("Applying optimal_dataset sample filtering")
    obj2 <- optimal_dataset(
      obj2,
      ident = ident_pref,
      sample = "sample",
      min_samples = processing_cfg$min_samples,
      max_samples = processing_cfg$max_samples,
      min_cells_sample = processing_cfg$min_cells_per_sample,
      verbose = TRUE
    )
  }

  if (is.null(obj2)) {
    warning("Skipping (insufficient samples): ", dataset_cfg$id, " / ", batch_val, " / ", condition_val)
    return(invisible(NULL))
  }

  dataset_key <- paste(
    sanitize_token(dataset_cfg$id),
    sanitize_token(paste(batch_vals[nzchar(batch_vals)], collapse = "-")),
    sanitize_token(paste(condition_vals[nzchar(condition_vals)], collapse = "-")),
    sep = "_"
  )
  dataset_key <- stringr::str_replace_all(dataset_key, "_+$", "")

  # ISC processed output
  ensure_dir(out_isc_dir)
  isc_path <- file.path(out_isc_dir, paste0(dataset_key, ".rds"))
  saveRDS(obj2, isc_path)

  # Sidecar metadata for downstream drivers
  meta <- list(
    dataset_key = dataset_key,
    dataset_id = dataset_cfg$id,
    dataset_reference = dataset_cfg$reference,
    batch = paste(batch_vals[nzchar(batch_vals)], collapse = ";"),
    condition = paste(condition_vals[nzchar(condition_vals)], collapse = ";"),
    idents = idents,
    isc = isTRUE(group_row$isc),
    batch_comparison = isTRUE(group_row$batch_comparison),
    perturbation_comparison = isTRUE(group_row$perturbation_comparison),
    label_transfer = isTRUE(group_row$label_transfer)
  )
  yaml::write_yaml(meta, paste0(isc_path, ".yaml"))

  message("Saved ISC processed dataset: ", isc_path)

  # Label-transfer query/reference split
  if (isTRUE(group_row$label_transfer)) {
    message_time("Preparing label-transfer query/reference split")

    counts <- SeuratObject::GetAssayData(obj2, assay = DefaultAssay(obj2), slot = "counts")

    mdlt <- obj2@meta.data
    mdlt$cell_type <- mdlt[[ident_pref]]
    mdlt <- mdlt %>%
      select(sample, cell_type) %>%
      mutate(dataset_id = dataset_key)

    split <- split_samples_query_reference(
      metadata = mdlt,
      sample_col = "sample",
      prop_query = label_transfer_cfg$prop_query_samples,
      seed = processing_cfg$seed
    )

    query_samples <- split$query_samples
    cells_query <- rownames(mdlt[mdlt$sample %in% query_samples, , drop = FALSE])
    cells_ref <- rownames(mdlt[!mdlt$sample %in% query_samples, , drop = FALSE])

    outdir <- file.path(out_lt_dir, dataset_key)
    ensure_dir(outdir)

    query <- list(counts = counts[, cells_query, drop = FALSE],
                  metadata = mdlt[cells_query, , drop = FALSE])
    reference <- list(counts = counts[, cells_ref, drop = FALSE],
                      metadata = mdlt[cells_ref, , drop = FALSE])

    saveRDS(query, file.path(outdir, "query.rds"))
    saveRDS(reference, file.path(outdir, "reference.rds"))

    message("Saved label-transfer split: ", outdir)
  }

  invisible(list(isc_path = isc_path, dataset_key = dataset_key))
}

main <- function() {
  stop_if_not_root()

  args <- parse_args(commandArgs(trailingOnly = TRUE))
  download <- tolower(args[["download"]] %||% "false") %in% c("true", "t", "1", "yes")

  specs <- read_specs("data_processing/specs_datasets.csv")
  registry <- read_dataset_registry("config/datasets.yaml")
  zenodo_doi <- read_zenodo_doi("config/datasets.yaml")
  registry_map <- registry_by_reference(registry)

  if (!is.na(zenodo_doi) && nzchar(zenodo_doi)) {
    message("Zenodo DOI: ", zenodo_doi)
  }

  params <- yaml::read_yaml("config/benchmark_parameters.yaml")
  seed <- params$seed %||% 22
  processing_cfg <- params$processing
  processing_cfg$seed <- seed
  label_transfer_cfg <- params$label_transfer

  zenodo_raw_dir <- "zenodo/raw"
  out_isc_dir <- "data/processed/isc"
  out_lt_dir <- "data/processed/label_transfer"

  # Expand + group specs to unique dataset-batch-condition entries.
  specs_groups <- specs %>%
    mutate(
      `Dataset reference` = normalize_text(`Dataset reference`),
      Condition = ifelse(is.na(Condition), "", Condition),
      Batch = ifelse(is.na(Batch), "", Batch),
      idents = strsplit(`# Annotation frameworks`, ",")
    ) %>%
    group_by(`Dataset reference`, Batch, Condition) %>%
    summarise(
      idents = list(unique(stringr::str_trim(unlist(idents)))),
      isc = any(`ISC-benchmarking` == "yes"),
      # Tasks 7–8 require combined objects (multiple batches/conditions). We therefore
      # do not flag the per-(Batch,Condition) objects for those tasks.
      batch_comparison = FALSE,
      perturbation_comparison = FALSE,
      label_transfer = any(`Label-Transfer Task` == "yes"),
      .groups = "drop"
    )

  # Additional combined datasets for manuscript Tasks 7–8.
  # - Task 7 (Batch comparison): combine *multiple batches* within a condition.
  # - Task 8 (Perturbation comparison): combine *multiple conditions* within a batch.
  specs_yes <- specs %>%
    mutate(
      `Dataset reference` = normalize_text(`Dataset reference`),
      Condition = ifelse(is.na(Condition), "", Condition),
      Batch = ifelse(is.na(Batch), "", Batch),
      idents = strsplit(`# Annotation frameworks`, ",")
    )

  batch_comp_groups <- specs_yes %>%
    filter(`Batch comparison` == "yes") %>%
    group_by(`Dataset reference`, Condition) %>%
    summarise(
      Batch = list(unique(Batch)),
      Condition = list(unique(Condition)),
      idents = list(unique(stringr::str_trim(unlist(idents)))),
      .groups = "drop"
    )

  pert_comp_groups <- specs_yes %>%
    filter(`Perturbation comparison` == "yes") %>%
    group_by(`Dataset reference`, Batch) %>%
    summarise(
      Batch = list(unique(Batch)),
      Condition = list(unique(Condition)),
      idents = list(unique(stringr::str_trim(unlist(idents)))),
      .groups = "drop"
    )

  # Download all configured files (if requested)
  if (isTRUE(download)) {
    message_time("Downloading Zenodo files (if missing)")
    for (ds in registry) {
      download_if_needed(ds, zenodo_raw_dir, download = TRUE)
    }
  }

  # Process each dataset-batch-condition group
  results <- list()
  for (i in seq_len(nrow(specs_groups))) {
    row <- specs_groups[i, ]
    ref <- normalize_text(row[["Dataset reference"]])

    if (!ref %in% names(registry_map)) {
      stop(
        "Dataset reference not found in config/datasets.yaml: ", ref,
        "\nTip: normalize whitespace in specs (or add a registry entry)."
      )
    }

    ds_cfg <- registry_map[[ref]]

    res <- process_one_spec_group(
      group_row = list(
        Batch = row$Batch,
        Condition = row$Condition,
        idents = row$idents[[1]],
        isc = row$isc,
        batch_comparison = row$batch_comparison,
        perturbation_comparison = row$perturbation_comparison,
        label_transfer = row$label_transfer
      ),
      dataset_cfg = ds_cfg,
      processing_cfg = processing_cfg,
      label_transfer_cfg = label_transfer_cfg,
      zenodo_raw_dir = zenodo_raw_dir,
      out_isc_dir = out_isc_dir,
      out_lt_dir = out_lt_dir,
      zenodo_doi = zenodo_doi
    )

    results[[i]] <- res
  }

  # Process combined datasets for Task 7 (batch comparison)
  for (i in seq_len(nrow(batch_comp_groups))) {
    row <- batch_comp_groups[i, ]
    ref <- normalize_text(row[["Dataset reference"]])
    batches <- row$Batch[[1]]
    conditions <- row$Condition[[1]]

    if (length(unique(batches[nzchar(batches)])) < 2) next
    if (!ref %in% names(registry_map)) {
      stop("Dataset reference not found in config/datasets.yaml: ", ref)
    }

    ds_cfg <- registry_map[[ref]]
    process_one_spec_group(
      group_row = list(
        Batch = batches,
        Condition = conditions,
        idents = row$idents[[1]],
        isc = FALSE,
        batch_comparison = TRUE,
        perturbation_comparison = FALSE,
        label_transfer = FALSE
      ),
      dataset_cfg = ds_cfg,
      processing_cfg = processing_cfg,
      label_transfer_cfg = label_transfer_cfg,
      zenodo_raw_dir = zenodo_raw_dir,
      out_isc_dir = out_isc_dir,
      out_lt_dir = out_lt_dir,
      zenodo_doi = zenodo_doi
    )
  }

  # Process combined datasets for Task 8 (perturbation comparison)
  for (i in seq_len(nrow(pert_comp_groups))) {
    row <- pert_comp_groups[i, ]
    ref <- normalize_text(row[["Dataset reference"]])
    batches <- row$Batch[[1]]
    conditions <- row$Condition[[1]]

    if (length(unique(conditions[nzchar(conditions)])) < 2) next
    if (!ref %in% names(registry_map)) {
      stop("Dataset reference not found in config/datasets.yaml: ", ref)
    }

    ds_cfg <- registry_map[[ref]]
    process_one_spec_group(
      group_row = list(
        Batch = batches,
        Condition = conditions,
        idents = row$idents[[1]],
        isc = FALSE,
        batch_comparison = FALSE,
        perturbation_comparison = TRUE,
        label_transfer = FALSE
      ),
      dataset_cfg = ds_cfg,
      processing_cfg = processing_cfg,
      label_transfer_cfg = label_transfer_cfg,
      zenodo_raw_dir = zenodo_raw_dir,
      out_isc_dir = out_isc_dir,
      out_lt_dir = out_lt_dir,
      zenodo_doi = zenodo_doi
    )
  }

  message_time("Dataset processing completed")
  invisible(results)
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

if (!interactive()) {
  main()
}
