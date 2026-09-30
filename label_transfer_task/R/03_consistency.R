# label_transfer_task/R/03_consistency.R --- Consistency + F1 metrics (scTypeEval)

source("R/00_utils.R")

# SCCAF is run the same way as in ISC_benchmark (external Python tool wrapper)
if (!exists("run_sample_agnostic_python_pipelines", mode = "function")) {
  source("../sample_agnostic_utils/run_python_pipelines.R")
}

# Cramer/Hotelling dataset-shift diagnostics (query vs reference)
if (!exists("run_scdiagnostics", mode = "function")) {
  source("../sample_agnostic_utils/scdiagnostics.R")
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(readr)
})

# Methods combined (as a product) into the single ISC summary score at
# aggregation time; kept out of compute_consistency_core() so referencing it
# there doesn't invalidate every already-cached per-branch consistency target.
LT_PRODUCT_CONSISTENCY_METHODS <- c(
  "silhouette | recip_classif:Match",
  "2label_silhouette | Pseudobulk:Cosine"
)

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

build_sc_for_scdiagnostics <- function(counts_matrix,
                                       metadata,
                                       ident,
                                       sample_col,
                                       min_samples) {
  metadata[[sample_col]] <- purge_label_local(metadata[[sample_col]])
  metadata[[ident]] <- purge_label_local(metadata[[ident]])

  sc <- scTypeEval::create_scTypeEval(
    matrix = counts_matrix,
    metadata = metadata,
    black_list = get_default_blacklist()
  )

  scTypeEval::run_processing_data(
    sc,
    ident = ident,
    sample = sample_col,
    min_samples = min_samples,
    verbose = FALSE
  )
}

compute_query_ref_scdiagnostics <- function(query_counts,
                                            query_metadata,
                                            ref_counts,
                                            ref_metadata,
                                            sample_col) {
  n_samples_query <- length(unique(query_metadata[[sample_col]]))
  n_samples_ref <- length(unique(ref_metadata[[sample_col]]))
  min_samples_query <- min(3, n_samples_query)
  min_samples_ref <- min(3, n_samples_ref)

  query_sc <- build_sc_for_scdiagnostics(
    counts_matrix = query_counts,
    metadata = query_metadata,
    ident = "pred_labels",
    sample_col = sample_col,
    min_samples = min_samples_query
  )

  ref_sc <- build_sc_for_scdiagnostics(
    counts_matrix = ref_counts,
    metadata = ref_metadata,
    ident = "true_labels",
    sample_col = sample_col,
    min_samples = min_samples_ref
  )

  run_scdiagnostics(query = query_sc, ref = ref_sc)
}

compute_sccaf_scores <- function(sc,
                                 ident,
                                 sample_col,
                                 min_samples,
                                 sccaf_n = 100,
                                 file_prefix = NULL) {
  sc_sccaf <- scTypeEval::run_processing_data(
    sc,
    ident = ident,
    aggregation = "single-cell",
    sample = sample_col,
    min_samples = min_samples,
    verbose = FALSE
  )

  if (is.null(file_prefix)) {
    file_prefix <- paste0("lt_sccaf_", format(Sys.time(), "%Y%m%d_%H%M%OS3"))
  }

  py_results <- run_sample_agnostic_python_pipelines(
    scTypeEval = sc_sccaf,
    pipelines = list(sccaf = pipeline_spec_sccaf(cluster_key = ident, n = sccaf_n)),
    tmp_dir = file.path(lt_consistency_dir(), "external_tmp"),
    file_prefix = file_prefix,
    continue_on_error = TRUE,
    cleanup = TRUE
  )

  py_results$sccaf
}

compute_consistency_core <- function(counts_matrix,
                                     metadata,
                                     ident,
                                     sample_col = "sample",
                                     method_diss = c("Pseudobulk:Cosine", "recip_classif:Match"),
                                     consistency_metric = c("silhouette",
                                                            "2label_silhouette",
                                                            "MetaNeighbor_Supervised",
                                                              "nsa_cLISI"
                                    ),
                                     cons_methods = c(
                                       "silhouette | recip_classif:Match",
                                       "2label_silhouette | Pseudobulk:Cosine",
                                       "MetaNeighbor_Supervised | NA",
                                       "nsa_cLISI | NA"
                                     ),
                                     ncores = 1,
                                     run_sccaf = TRUE,
                                     sccaf_n = 100,
                                     file_prefix = NULL) {
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

  n_samples <- length(unique(metadata[[sample_col]]))
  if (n_samples < 2) {
    stop("Need >= 2 samples to compute consistency; found ", n_samples)
  }
  min_samples <- min(3, n_samples)

  sc <- scTypeEval::run_processing_data(
    sc,
    ident = ident,
    sample = sample_col,
    min_samples = min_samples,
    verbose = FALSE
  )

  sc <- scTypeEval::run_hvg(
    sc,
    ngenes = 2000,
    aggregation = "single-cell",
    ncores = ncores,
    verbose = FALSE
  )

  sc <- scTypeEval::run_pca(sc, verbose = FALSE)

  for (mdiss in method_diss) {
    sc <- scTypeEval::run_dissimilarity(
      sc,
      method = mdiss,
      ncores = ncores,
      verbose = FALSE
    )
  }

  cons <- scTypeEval::get_consistency(sc,
                                      consistency_metric = consistency_metric,
                                      verbose = FALSE) %>%
    dplyr::rename(cell_type = celltype) %>%
    dplyr::mutate(method_type = paste(consistency_metric, dissimilarity_method, sep = " | ")) %>%
    dplyr::filter(method_type %in% cons_methods) %>%
    dplyr::select(-consistency_metric, -dissimilarity_method) %>%
    tidyr::pivot_wider(names_from = method_type, values_from = measure)

  if (isTRUE(run_sccaf)) {
    sccaf_scores <- tryCatch(
      compute_sccaf_scores(
        sc = sc,
        ident = ident,
        sample_col = sample_col,
        min_samples = min_samples,
        sccaf_n = sccaf_n,
        file_prefix = file_prefix
      ),
      error = function(e) {
        message("[SCCAF] Skipped: ", conditionMessage(e))
        NULL
      }
    )

    cons$SCCAF <- NA_real_
    if (!is.null(sccaf_scores) && nrow(sccaf_scores) > 0) {
      sccaf_scores <- sccaf_scores %>%
        dplyr::transmute(cell_type = as.character(celltype), SCCAF = as.numeric(score))
      cons <- cons %>%
        dplyr::select(-SCCAF) %>%
        dplyr::left_join(sccaf_scores, by = "cell_type")
    }
  }

  cons
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
                                         ncores = 1,
                                         reference_dataset_id = NULL,
                                         query_dataset_id = NULL) {
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

  cons <- tryCatch(
    compute_consistency_core(
      counts_matrix = counts,
      metadata = md,
      ident = "pred_labels",
      sample_col = sample_col,
      ncores = ncores,
      file_prefix = sprintf("%s_%s_rep%d_query", dataset_id, classifier_name, rep)
    ),
    error = function(e) {
      message(sprintf("[SKIP] Query consistency failed for '%s' / '%s' (rep %d): %s",
                      dataset_id, classifier_name, rep, conditionMessage(e)))
      return(invisible(NULL))
    }
  )
  if (is.null(cons)) return(invisible(NULL))

  cons <- add_f1_one_vs_rest(cons, pred_labels = md$pred_labels, true_labels = md$true_labels)

  # Cramer/Hotelling: do query cells assigned each predicted label distributionally
  # match reference cells truly labeled with that same cell type?
  ref_path <- file.path(dataset_dir, "reference.rds")
  if (file.exists(ref_path)) {
    scdiag <- tryCatch({
      ref <- readRDS(ref_path)
      ref_counts <- ref$counts
      ref_cell_ids <- colnames(ref_counts)
      ref_md_raw <- as.data.frame(ref$metadata)

      if (is.null(rownames(ref_md_raw)) ||
          all(rownames(ref_md_raw) == as.character(seq_len(nrow(ref_md_raw))))) {
        if (nrow(ref_md_raw) != length(ref_cell_ids)) {
          stop("Reference metadata/counts mismatch for ", dataset_id, ": ",
               nrow(ref_md_raw), " rows vs ", length(ref_cell_ids), " cells")
        }
        rownames(ref_md_raw) <- ref_cell_ids
        ref_md <- ref_md_raw
      } else {
        ref_md <- ref_md_raw[ref_cell_ids, , drop = FALSE]
      }
      ref_md$true_labels <- ref_md$cell_type

      compute_query_ref_scdiagnostics(
        query_counts = counts,
        query_metadata = md,
        ref_counts = ref_counts,
        ref_metadata = ref_md,
        sample_col = sample_col
      )
    }, error = function(e) {
      message(sprintf("[SKIP] Cramer/Hotelling failed for '%s' / '%s' (rep %d): %s",
                      dataset_id, classifier_name, rep, conditionMessage(e)))
      return(invisible(NULL))
    })

    if (!is.null(scdiag) && nrow(scdiag) > 0) {
      scdiag <- scdiag %>%
        dplyr::transmute(cell_type = as.character(celltype), cramer, hotelling)
      cons <- cons %>%
        dplyr::left_join(scdiag, by = "cell_type")
    }
  }

  cons <- cons %>%
    dplyr::mutate(
      dataset_id = dataset_id,
      classifier = classifier_name,
      replicate = rep,
      split = "query"
    )

  if (!is.null(reference_dataset_id)) {
    cons$reference_dataset_id <- reference_dataset_id
  }
  if (!is.null(query_dataset_id)) {
    cons$query_dataset_id <- query_dataset_id
  }

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
                                             ncores = 1,
                                             reference_dataset_id = NULL,
                                             query_dataset_id = NULL) {
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

  cons <- tryCatch(
    compute_consistency_core(
      counts_matrix = counts,
      metadata = md,
      ident = "true_labels",
      sample_col = sample_col,
      ncores = ncores,
      file_prefix = sprintf("%s_rep%d_reference", dataset_id, rep)
    ),
    error = function(e) {
      message(sprintf("[SKIP] Reference consistency failed for '%s' (rep %d): %s",
                      dataset_id, rep, conditionMessage(e)))
      return(invisible(NULL))
    }
  )
  if (is.null(cons)) return(invisible(NULL))

  cons <- cons %>%
    dplyr::mutate(
      dataset_id = dataset_id,
      classifier = "ground_truth",
      replicate = rep,
      split = "reference"
    )

  if (!is.null(reference_dataset_id)) {
    cons$reference_dataset_id <- reference_dataset_id
  }
  if (!is.null(query_dataset_id)) {
    cons$query_dataset_id <- query_dataset_id
  }

  out_file <- file.path(
    output_dir,
    sprintf("%s_rep%d_reference_ground_truth_consistency.rds", dataset_id, rep)
  )
  saveRDS(cons, out_file)

  invisible(out_file)
}

# Between-dataset pairs sharing the same reference_dataset_id have identical
# reference.rds contents, so compute reference consistency once per unique
# (reference_dataset_id, rep) and fan the result out to every pair below,
# instead of recomputing per pair_id.
compute_lt_between_unique_reference_consistency <- function(reference_dataset_id,
                                                             rep,
                                                             pairs,
                                                             data_dir,
                                                             sample_col = "sample",
                                                             ncores = 1) {
  representative_pair_id <- pairs$pair_id[pairs$reference_dataset_id == reference_dataset_id][1]
  if (is.na(representative_pair_id)) return(invisible(NULL))

  ref_path <- file.path(data_dir, representative_pair_id, "reference.rds")
  if (!file.exists(ref_path)) return(invisible(NULL))

  ref <- readRDS(ref_path)
  counts <- ref$counts
  cell_ids <- colnames(counts)
  md_raw <- as.data.frame(ref$metadata)

  if (is.null(rownames(md_raw)) || all(rownames(md_raw) == as.character(seq_len(nrow(md_raw))))) {
    if (nrow(md_raw) != length(cell_ids)) {
      stop("Reference metadata/counts mismatch for ", reference_dataset_id, ": ",
           nrow(md_raw), " rows vs ", length(cell_ids), " cells")
    }
    rownames(md_raw) <- cell_ids
    md <- md_raw
  } else {
    md <- md_raw[cell_ids, , drop = FALSE]
  }
  md$true_labels <- md$cell_type

  cons <- tryCatch(
    compute_consistency_core(
      counts_matrix = counts,
      metadata = md,
      ident = "true_labels",
      sample_col = sample_col,
      ncores = ncores,
      file_prefix = sprintf("%s_rep%d_reference_shared", reference_dataset_id, rep)
    ),
    error = function(e) {
      message(sprintf("[SKIP] Shared reference consistency failed for '%s' (rep %d): %s",
                      reference_dataset_id, rep, conditionMessage(e)))
      return(invisible(NULL))
    }
  )
  if (is.null(cons)) return(invisible(NULL))

  cons$reference_dataset_id <- reference_dataset_id
  cons$replicate <- rep
  cons
}

# Relabels the shared reference-consistency table (dataset_id/query_dataset_id)
# per pair_id and writes one file per pair, so downstream joins on dataset_id
# (== pair_id) keep working while the expensive computation ran only once.
write_lt_between_reference_consistency_outputs <- function(unique_cons,
                                                            unique_grid,
                                                            pairs,
                                                            output_dir = NULL) {
  if (is.null(output_dir)) output_dir <- lt_between_consistency_dir()
  ensure_dir(output_dir)

  out_files <- character(0)
  for (i in seq_len(nrow(unique_grid))) {
    cons_base <- unique_cons[[i]]
    if (is.null(cons_base)) next

    ref_id <- unique_grid$reference_dataset_id[i]
    rep <- unique_grid$replicate[i]

    matched_pairs <- pairs[pairs$reference_dataset_id == ref_id, , drop = FALSE]
    for (j in seq_len(nrow(matched_pairs))) {
      pair_id <- matched_pairs$pair_id[j]
      query_dataset_id <- matched_pairs$query_dataset_id[j]

      cons <- cons_base %>%
        dplyr::mutate(
          dataset_id = pair_id,
          classifier = "ground_truth",
          replicate = rep,
          split = "reference",
          query_dataset_id = query_dataset_id
        )

      out_file <- file.path(
        output_dir,
        sprintf("%s_rep%d_reference_ground_truth_consistency.rds", pair_id, rep)
      )
      saveRDS(cons, out_file)
      out_files <- c(out_files, out_file)
    }
  }

  out_files
}

aggregate_lt_consistency_results <- function(consistency_dir, output_file) {
  ensure_dir(dirname(output_file))
  files <- list.files(consistency_dir, pattern = "\\.rds$", full.names = TRUE)
  if (length(files) == 0) {
    warning("No consistency result files found in ", consistency_dir)
    return(invisible(NULL))
  }

  combined <- purrr::map_df(files, readRDS)

  # "product" (combined ISC score) was never written per-branch, so derive it
  # here from the core method columns instead of recomputing every branch.
  metric_cols <- intersect(LT_PRODUCT_CONSISTENCY_METHODS, colnames(combined))
  combined$product <- vapply(seq_len(nrow(combined)), function(i) {
    vals <- as.numeric(combined[i, metric_cols, drop = TRUE])
    vals <- vals[!is.na(vals)]
    if (length(vals) == 0) return(NA_real_)
    prod(vals)
  }, numeric(1))

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
