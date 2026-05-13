# Meta_study/R/00_utils.R --- Utility helpers for the Meta-study targets pipeline

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

resolve_meta_path <- function(path_value, must_work = TRUE) {
  if (is.null(path_value) || !nzchar(path_value)) {
    return(path_value)
  }

  if (startsWith(path_value, "/")) {
    return(normalizePath(path_value, mustWork = must_work))
  }

  normalizePath(proj_path(path_value), mustWork = must_work)
}

load_meta_params <- function(config_path = "config/meta_study_parameters.yaml") {
  params <- yaml::read_yaml(config_path)

  params$paths$data_local_dir <- resolve_meta_path(params$paths$data_local_dir, must_work = FALSE)
  params$paths$data_nas_root <- resolve_meta_path(params$paths$data_nas_root, must_work = FALSE)
  params$paths$output_scTypeEval_dir <- resolve_meta_path(params$paths$output_scTypeEval_dir, must_work = FALSE)
  params$paths$scTypeEval_source_dir <- resolve_meta_path(params$paths$scTypeEval_source_dir, must_work = FALSE)
  params$paths$default_black_list_rdata <- resolve_meta_path(params$paths$default_black_list_rdata, must_work = TRUE)

  ensure_dir(params$paths$output_scTypeEval_dir)
  params
}

load_meta_datasets <- function(config_path = "config/datasets_metadata.yaml") {
  datasets <- yaml::read_yaml(config_path)$datasets
  if (is.null(datasets) || length(datasets) == 0) {
    stop("No datasets found in: ", config_path)
  }

  tibble::tibble(
    name = vapply(datasets, function(x) x$name %||% NA_character_, character(1)),
    ds_name = vapply(datasets, function(x) x$ds_name %||% NA_character_, character(1)),
    sample_col = vapply(datasets, function(x) x$sample_col %||% NA_character_, character(1)),
    low_res_ct_col = vapply(datasets, function(x) x$low_res_ct_col %||% NA_character_, character(1)),
    hi_res_ct_col = vapply(datasets, function(x) x$hi_res_ct_col %||% NA_character_, character(1)),
    nas_rel_path = vapply(datasets, function(x) x$nas_rel_path %||% NA_character_, character(1))
  )
}

get_meta_dissimilarity_methods <- function(params) {
  methods <- params$dissimilarity_methods %||% character(0)
  methods <- as.character(unlist(methods, use.names = FALSE))
  methods <- methods[nzchar(methods)]
  if (length(methods) == 0) {
    stop("No dissimilarity methods configured in meta_study_parameters.yaml")
  }
  methods
}

ensure_scTypeEval_api <- function(params) {
  required_symbols <- c(
    "load_singleCell_object",
    "create.scTypeEval",
    "Run.ProcessingData",
    "Run.HVG",
    "Run.PCA",
    "Run.Dissimilarity"
  )

  if (requireNamespace("scTypeEval", quietly = TRUE)) {
    suppressPackageStartupMessages(library(scTypeEval))
  }

  missing_symbols <- required_symbols[!vapply(required_symbols, exists, logical(1), mode = "function")]

  if (length(missing_symbols) > 0) {
    src_dir <- params$paths$scTypeEval_source_dir
    r_scripts <- list.files(src_dir, pattern = "\\.R$", full.names = TRUE)
    if (length(r_scripts) == 0) {
      stop("No R scripts found in scTypeEval source dir: ", src_dir)
    }
    for (script in r_scripts) {
      source(script)
    }
  }

  missing_symbols <- required_symbols[!vapply(required_symbols, exists, logical(1), mode = "function")]
  if (length(missing_symbols) > 0) {
    stop("Missing required scTypeEval functions after loading: ", paste(missing_symbols, collapse = ", "))
  }

  invisible(TRUE)
}

load_black_list_vector <- function(params) {
  bl_env <- new.env(parent = emptyenv())
  load(params$paths$default_black_list_rdata, envir = bl_env)

  if (!exists("black.list", envir = bl_env, inherits = FALSE)) {
    stop("Object 'black.list' was not found in file: ", params$paths$default_black_list_rdata)
  }

  black.list <- get("black.list", envir = bl_env, inherits = FALSE)
  unlist(list(black.list$TCR, black.list$Immunoglobulins, black.list$Ygenes), use.names = FALSE)
}

build_meta_jobs <- function(datasets_tbl) {
  jobs <- list()

  for (i in seq_len(nrow(datasets_tbl))) {
    row <- datasets_tbl[i, , drop = FALSE]
    ident_cols <- c(row$low_res_ct_col, row$hi_res_ct_col)
    ident_cols <- ident_cols[!is.na(ident_cols) & nzchar(ident_cols)]

    if (length(ident_cols) == 0) {
      next
    }

    for (ident_col in ident_cols) {
      jobs[[length(jobs) + 1]] <- tibble::tibble(
        name = row$name,
        ds_name = row$ds_name,
        sample_col = row$sample_col,
        ident_col = ident_col,
        nas_rel_path = row$nas_rel_path
      )
    }
  }

  if (length(jobs) == 0) {
    stop("No dataset/annotation jobs were generated from config/datasets_metadata.yaml")
  }

  dplyr::bind_rows(jobs)
}

resolve_dataset_file <- function(ds_row, params) {
  local_file <- file.path(params$paths$data_local_dir, paste0(ds_row$ds_name, ".rds"))
  if (file.exists(local_file)) {
    return(local_file)
  }

  nas_file <- file.path(params$paths$data_nas_root, ds_row$nas_rel_path)
  if (file.exists(nas_file)) {
    return(nas_file)
  }

  stop(
    "Dataset file not found for ", ds_row$ds_name,
    ". Checked local: ", local_file,
    " and NAS: ", nas_file
  )
}

run_meta_scTypeeval_job <- function(ds_row, ident_col, params, black_list, dissimilarity_methods) {
  output_file <- file.path(
    params$paths$output_scTypeEval_dir,
    paste0(ds_row$ds_name, "__", ident_col, ".rds")
  )

  if (isTRUE(params$execution$skip_existing) && file.exists(output_file)) {
    message("[SKIP] ", basename(output_file), " already exists")
    return(output_file)
  }

  dataset_file <- resolve_dataset_file(ds_row, params)
  message("[RUN] ", ds_row$ds_name, " (", ident_col, ")")

  object <- load_singleCell_object(dataset_file)
  metadata <- object[["metadata"]]

  if (!(ds_row$sample_col %in% colnames(metadata))) {
    stop("Sample column not found in metadata: ", ds_row$sample_col)
  }
  if (!(ident_col %in% colnames(metadata))) {
    stop("Ident column not found in metadata: ", ident_col)
  }

  sub_meta <- data.frame(
    sample = metadata[[ds_row$sample_col]],
    stringsAsFactors = FALSE
  )
  sub_meta[[ident_col]] <- metadata[[ident_col]]
  rownames(sub_meta) <- rownames(metadata)
  sub_meta <- sub_meta[!is.na(sub_meta$sample) & !is.na(sub_meta[[ident_col]]), , drop = FALSE]

  if (nrow(sub_meta) == 0) {
    stop("No cells left after filtering NA sample/ident for ", ds_row$ds_name, " and ", ident_col)
  }

  sub_count_matrix <- object[["counts"]][, rownames(sub_meta), drop = FALSE]

  sc <- create.scTypeEval(
    matrix = sub_count_matrix,
    metadata = sub_meta,
    black.list = black_list
  )

  sc <- Run.ProcessingData(
    sc,
    ident = ident_col,
    sample = "sample",
    verbose = FALSE
  )

  sc <- Run.HVG(
    sc,
    ngenes = params$execution$hvg_genes,
    ncores = params$execution$ncores,
    verbose = FALSE
  )

  sc@data$`single-cell` <- NULL
  sc@counts <- as(Matrix::Matrix(0, nrow = 0, ncol = 0, sparse = TRUE), "dgCMatrix")

  sc <- Run.PCA(sc, verbose = FALSE)

  for (mdiss in dissimilarity_methods) {
    sc <- Run.Dissimilarity(
      sc,
      method = mdiss,
      ncores = params$execution$ncores,
      verbose = FALSE
    )
  }

  saveRDS(sc, output_file)

  rm(object, metadata, sub_meta, sub_count_matrix, sc)
  gc()

  output_file
}

write_meta_summary <- function(output_files, params) {
  output_files <- unlist(output_files, use.names = FALSE)
  output_files <- output_files[file.exists(output_files)]

  summary_df <- tibble::tibble(file = basename(output_files))
  summary_df$dataset_id <- sub("__.*$", "", tools::file_path_sans_ext(summary_df$file))
  summary_df$annotation <- sub("^.*__", "", tools::file_path_sans_ext(summary_df$file))

  summary_file <- file.path(params$paths$output_scTypeEval_dir, "meta_study_scTypeEval_summary.csv")
  utils::write.csv(summary_df, summary_file, row.names = FALSE)
  summary_file
}
