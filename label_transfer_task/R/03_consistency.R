# label_transfer_task/R/03_consistency.R --- Consistency + F1 metrics (scTypeEval)

source("R/00_utils.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(readr)
})

purge_label_local <- function(label) {
  # Mirror scTypeEval::purge_label() behavior without relying on :::
  label <- as.character(label)
  label <- gsub(" |_|[+]|-", ".", label)
  label <- gsub(",", "", label)
  label
}

get_default_blacklist <- function() {
  if (!requireNamespace("scTypeEval", quietly = TRUE)) {
    stop("Package 'scTypeEval' is required for consistency computation")
  }

  # Provided by scTypeEval data
  data("black_list", package = "scTypeEval", envir = environment())

  unlist(list(
    black_list$TCR,
    black_list$Immunoglobulins,
    black_list$Ygenes
  ))
}

compute_consistency_core <- function(counts_matrix,
                                     metadata,
                                     ident,
                                     sample_col = "sample",
                                     method_diss = c("Pseudobulk:Cosine", "recip_classif:Match"),
                                     cons_methods = c(
                                       "silhouette | recip_classif:Match",
                                       "2label_silhouette | Pseudobulk:Cosine"
                                     ),
                                     ncores = 1) {
  if (!requireNamespace("scTypeEval", quietly = TRUE)) {
    stop("Package 'scTypeEval' is required for consistency computation")
  }

  metadata[[sample_col]] <- purge_label_local(metadata[[sample_col]])
  metadata[[ident]] <- purge_label_local(metadata[[ident]])

  sc <- scTypeEval::create_scTypeEval(
    matrix = counts_matrix,
    metadata = metadata,
    black_list = get_default_blacklist()
  )

  hvg_ident <- if ("true_labels" %in% colnames(metadata)) "true_labels" else ident

  n_samples <- length(unique(metadata[[sample_col]]))
  if (n_samples < 2) {
    stop("Need >= 2 samples to compute consistency; found ", n_samples)
  }
  min_samples <- min(3, n_samples)

  sc_tmp <- scTypeEval::run_processing_data(
    sc,
    ident = hvg_ident,
    aggregation = "single-cell",
    sample = sample_col,
    min_samples = min_samples,
    verbose = FALSE
  )

  sc_tmp <- scTypeEval::run_hvg(
    sc_tmp,
    ngenes = 2000,
    ncores = ncores,
    verbose = FALSE
  )

  hvg <- sc_tmp@gene_lists

  sc_proc <- scTypeEval::run_processing_data(
    sc,
    ident = ident,
    aggregation = "pseudobulk",
    sample = sample_col,
    min_samples = min_samples,
    verbose = FALSE
  )

  sc_proc <- scTypeEval::add_gene_list(sc_proc, gene_list = hvg)
  sc_proc <- scTypeEval::run_pca(sc_proc, verbose = FALSE)

  for (mdiss in method_diss) {
    sc_proc <- scTypeEval::run_dissimilarity(
      sc_proc,
      method = mdiss,
      ncores = ncores,
      verbose = FALSE
    )
  }

  scTypeEval::get_consistency(sc_proc) %>%
    dplyr::rename(cell_type = celltype) %>%
    dplyr::mutate(method_type = paste(consistency_metric, dissimilarity_method, sep = " | ")) %>%
    dplyr::filter(method_type %in% cons_methods) %>%
    dplyr::select(-consistency_metric, -dissimilarity_method) %>%
    tidyr::pivot_wider(names_from = method_type, values_from = measure) %>%
    dplyr::mutate(product = .data[[cons_methods[1]]] * .data[[cons_methods[2]]])
}

add_f1_one_vs_rest <- function(cons_table, pred_labels, true_labels) {
  pred_vector <- purge_label_local(pred_labels)
  true_vector <- purge_label_local(true_labels)

  cell_types <- unique(as.character(true_vector))
  cell_types <- cell_types[!is.na(cell_types)]

  per_type_metrics <- lapply(cell_types, function(ct) {
    pred_binary <- as.character(pred_vector) == ct
    true_binary <- as.character(true_vector) == ct

    valid <- !(is.na(pred_binary) | is.na(true_binary))
    pred_binary <- pred_binary[valid]
    true_binary <- true_binary[valid]

    tp <- sum(pred_binary & true_binary)
    fp <- sum(pred_binary & !true_binary)
    fn <- sum(!pred_binary & true_binary)
    tn <- sum(!pred_binary & !true_binary)

    acc <- (tp + tn) / (tp + fp + fn + tn)
    f1_val <- if ((tp + fp + fn) > 0) 2 * tp / (2 * tp + fp + fn) else 0

    data.frame(cell_type = ct, accuracy = acc, f1 = f1_val)
  })

  per_type_df <- do.call(rbind, per_type_metrics)
  rownames(per_type_df) <- NULL

  cons_table %>%
    dplyr::left_join(per_type_df, by = "cell_type")
}

compute_lt_query_consistency <- function(dataset_id,
                                         classifier_name,
                                         rep,
                                         result_path = NULL,
                                         data_dir = NULL,
                                         results_dir = NULL,
                                         output_dir = NULL,
                                         sample_col = "sample",
                                         ncores = 1) {
  if (is.null(data_dir)) data_dir <- lt_data_processed_dir(rep)
  if (is.null(results_dir)) results_dir <- lt_raw_results_dir()
  if (is.null(output_dir)) output_dir <- lt_consistency_dir()

  dataset_dir <- file.path(data_dir, dataset_id)
  query_path <- file.path(dataset_dir, "query.rds")

  if (is.null(result_path)) {
    result_path <- file.path(
      results_dir,
      sprintf("%s_%s_rep%d.rds", dataset_id, classifier_name, rep)
    )
  }

  if (!file.exists(query_path) || !file.exists(result_path)) {
    return(invisible(NULL))
  }

  query <- readRDS(query_path)
  counts <- query$counts
  run_df <- readRDS(result_path)

  cell_ids <- colnames(counts)

  if ("cell_id" %in% colnames(run_df)) {
    pred <- run_df$prediction
    names(pred) <- run_df$cell_id
    pred <- pred[cell_ids]
  } else {
    pred <- run_df$prediction
    if (length(pred) != length(cell_ids)) {
      stop("Prediction length mismatch for ", dataset_id, ": ",
           length(pred), " vs ", length(cell_ids))
    }
    names(pred) <- cell_ids
  }

  md_raw <- as.data.frame(query$metadata)

  if (is.null(rownames(md_raw)) || all(rownames(md_raw) == as.character(seq_len(nrow(md_raw))))) {
    if (nrow(md_raw) != length(cell_ids)) {
      stop("Metadata/Counts mismatch for ", dataset_id, ": ",
           nrow(md_raw), " rows vs ", length(cell_ids), " cells")
    }
    rownames(md_raw) <- cell_ids
    md <- md_raw
  } else {
    md <- md_raw[cell_ids, , drop = FALSE]
  }

  md$pred_labels <- unname(pred[cell_ids])
  md$true_labels <- md$cell_type

  cons <- compute_consistency_core(
    counts_matrix = counts,
    metadata = md,
    ident = "pred_labels",
    sample_col = sample_col,
    ncores = ncores
  )

  cons <- add_f1_one_vs_rest(cons, pred_labels = md$pred_labels, true_labels = md$true_labels)

  cons <- cons %>%
    dplyr::mutate(
      dataset_id = dataset_id,
      classifier = classifier_name,
      replicate = rep,
      split = "query"
    )

  out_file <- file.path(
    output_dir,
    sprintf("%s_%s_rep%d_query_consistency.rds", dataset_id, classifier_name, rep)
  )
  saveRDS(cons, out_file)

  invisible(out_file)
}

compute_lt_reference_consistency <- function(dataset_id,
                                             rep,
                                             data_dir = NULL,
                                             output_dir = NULL,
                                             sample_col = "sample",
                                             ncores = 1) {
  if (is.null(data_dir)) data_dir <- lt_data_processed_dir(rep)
  if (is.null(output_dir)) output_dir <- lt_consistency_dir()

  dataset_dir <- file.path(data_dir, dataset_id)
  ref_path <- file.path(dataset_dir, "reference.rds")

  if (!file.exists(ref_path)) {
    return(invisible(NULL))
  }

  ref <- readRDS(ref_path)
  counts <- ref$counts
  cell_ids <- colnames(counts)

  md_raw <- as.data.frame(ref$metadata)

  if (is.null(rownames(md_raw)) || all(rownames(md_raw) == as.character(seq_len(nrow(md_raw))))) {
    if (nrow(md_raw) != length(cell_ids)) {
      stop("Reference metadata/counts mismatch for ", dataset_id, ": ",
           nrow(md_raw), " rows vs ", length(cell_ids), " cells")
    }
    rownames(md_raw) <- cell_ids
    md <- md_raw
  } else {
    md <- md_raw[cell_ids, , drop = FALSE]
  }

  md$true_labels <- md$cell_type

  cons <- compute_consistency_core(
    counts_matrix = counts,
    metadata = md,
    ident = "true_labels",
    sample_col = sample_col,
    ncores = ncores
  )

  cons <- cons %>%
    dplyr::mutate(
      dataset_id = dataset_id,
      classifier = "ground_truth",
      replicate = rep,
      split = "reference"
    )

  out_file <- file.path(
    output_dir,
    sprintf("%s_rep%d_reference_ground_truth_consistency.rds", dataset_id, rep)
  )
  saveRDS(cons, out_file)

  invisible(out_file)
}

aggregate_lt_consistency_results <- function(consistency_dir, output_file) {
  ensure_dir(dirname(output_file))

  files <- list.files(consistency_dir, pattern = "\\.rds$", full.names = TRUE)
  if (length(files) == 0) {
    warning("No consistency result files found in ", consistency_dir)
    return(invisible(NULL))
  }

  combined <- purrr::map_df(files, readRDS)

  summary <- combined %>%
    dplyr::group_by(dataset_id, classifier, replicate, split) %>%
    dplyr::summarise(
      mean_product = mean(product, na.rm = TRUE),
      macro_f1 = {
        m <- mean(f1, na.rm = TRUE)
        if (is.nan(m)) NA_real_ else m
      },
      .groups = "drop"
    )

  detailed_file <- sub("\\.csv$", "_detailed.csv", output_file)
  summary_file <- sub("\\.csv$", "_summary.csv", output_file)

  readr::write_csv(combined, detailed_file)
  readr::write_csv(summary, summary_file)

  invisible(summary_file)
}
