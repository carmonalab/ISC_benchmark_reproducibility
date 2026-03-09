#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(yaml)
  library(dplyr)
})

source("utils/project_paths.R")
source("utils/cli_utils.R")

stop_if_not_root <- function() {
  if (!dir.exists("data/processed/label_transfer") || !file.exists("label_transfer_task/classifiers/classifiers.R")) {
    stop("Run from repository root; missing processed label_transfer data or classifier code.")
  }
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

normalize_text <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x <- gsub("[\u00A0\u200B\u200C\u200D\uFEFF]", " ", x, perl = TRUE)
  x <- gsub("[[:space:]]+", " ", x, perl = TRUE)
  trimws(x)
}

read_zenodo_doi <- function(path = "config/datasets.yaml") {
  if (!file.exists(path)) return(NA_character_)
  cfg <- yaml::read_yaml(path)
  doi <- normalize_text(cfg$zenodo$doi)
  if (!nzchar(doi) || identical(doi, "TODO")) return(NA_character_)
  doi
}

accuracy_metrics <- function(truth, estimate) {
  truth <- as.character(truth)
  estimate <- as.character(estimate)

  keep <- !is.na(truth) & !is.na(estimate)
  truth <- truth[keep]
  estimate <- estimate[keep]

  if (length(truth) == 0) {
    return(data.frame(accuracy = NA_real_, balanced_accuracy = NA_real_))
  }

  acc <- mean(truth == estimate)

  # Balanced accuracy: mean per-class recall
  classes <- sort(unique(truth))
  recalls <- vapply(classes, function(cls) {
    idx <- truth == cls
    if (sum(idx) == 0) return(NA_real_)
    mean(estimate[idx] == cls)
  }, numeric(1))

  data.frame(
    accuracy = acc,
    balanced_accuracy = mean(recalls, na.rm = TRUE)
  )
}

classify_one_dataset <- function(dataset_dir) {
  dataset_id <- basename(dataset_dir)
  query_file <- file.path(dataset_dir, "query.rds")
  reference_file <- file.path(dataset_dir, "reference.rds")

  if (!file.exists(query_file) || !file.exists(reference_file)) {
    warning("Missing query/reference in ", dataset_dir)
    return(NULL)
  }

  query <- readRDS(query_file)
  reference <- readRDS(reference_file)

  ref_counts <- reference$counts
  ref_labels <- reference$metadata$cell_type
  query_counts <- query$counts
  query_md <- query$metadata

  if (is.null(query_md$cell_type)) {
    stop("Query metadata must include cell_type for evaluation: ", dataset_id)
  }

  # Load classifier functions
  source("label_transfer_task/classifiers/classifiers.R")

  classifiers <- list(
    SingleR = classify_SingleR,
    LogisticRegression = classify_LogisticRegression,
    XGBoost = classify_XGBoost,
    MLP = classify_MLP,
    RandomForest = classify_RandomForest,
    SVM = classify_SVM,
    kNN = classify_kNN,
    NaiveBayes = classify_NaiveBayes,
    LDA = classify_LDA,
    SeuratTransfer = classify_SeuratTransfer,
    DecisionTree = classify_DecisionTree,
    scPred = classify_scPred,
    Random = classify_Random
  )

  predictions <- list()
  metrics <- list()

  for (clf in names(classifiers)) {
    message("[", dataset_id, "] ", clf)
    pred <- tryCatch(
      classifiers[[clf]](ref_counts, ref_labels, query_counts),
      error = function(e) {
        warning("Classifier failed: ", clf, " (", dataset_id, "): ", e$message)
        NA
      }
    )

    if (length(pred) != ncol(query_counts)) {
      pred <- rep(NA, ncol(query_counts))
    }

    predictions[[clf]] <- pred
    m <- accuracy_metrics(query_md$cell_type, pred)
    metrics[[clf]] <- cbind(
      dataset_id = dataset_id,
      classifier = clf,
      n_query_cells = ncol(query_counts),
      n_ref_cells = ncol(ref_counts),
      n_cell_types = dplyr::n_distinct(ref_labels),
      m
    )
  }

  # Ensemble voting (exclude Random)
  message("[", dataset_id, "] Ensemble")
  preds_for_ensemble <- predictions[names(predictions) != "Random"]
  ensemble_pred <- tryCatch(
    classify_Ensemble(preds_for_ensemble),
    error = function(e) {
      warning("Ensemble failed: ", dataset_id, ": ", e$message)
      rep(NA, ncol(query_counts))
    }
  )
  predictions$Ensemble <- ensemble_pred
  m <- accuracy_metrics(query_md$cell_type, ensemble_pred)
  metrics$Ensemble <- cbind(
    dataset_id = dataset_id,
    classifier = "Ensemble",
    n_query_cells = ncol(query_counts),
    n_ref_cells = ncol(ref_counts),
    n_cell_types = dplyr::n_distinct(ref_labels),
    m
  )

  # Save predictions as updated query metadata
  updated_md <- query_md
  for (clf in names(predictions)) {
    updated_md[[paste0("pred_", clf)]] <- predictions[[clf]]
  }

  saveRDS(updated_md, file.path(dataset_dir, "predictions.rds"))

  dplyr::bind_rows(metrics)
}

main <- function() {
  stop_if_not_root()

  args <- parse_args(commandArgs(trailingOnly = TRUE))
  out_dir <- args[["out"]] %||% "results/label_transfer"
  ensure_dir(out_dir)

  zenodo_doi <- read_zenodo_doi("config/datasets.yaml")
  if (!is.na(zenodo_doi) && nzchar(zenodo_doi)) {
    message("Zenodo DOI: ", zenodo_doi)
    yaml::write_yaml(
      list(
        zenodo_doi = zenodo_doi,
        script = "label_transfer_task/run_label_transfer_benchmark.R"
      ),
      file.path(out_dir, "run_metadata.yaml")
    )
  }

  dataset_dirs <- list.dirs("data/processed/label_transfer", recursive = FALSE, full.names = TRUE)
  dataset_dirs <- dataset_dirs[basename(dataset_dirs) != ".DS_Store"]
  if (length(dataset_dirs) == 0) {
    stop("No label-transfer datasets found under data/processed/label_transfer/. Run processing first.")
  }

  message_time(paste("Running label transfer benchmark on", length(dataset_dirs), "datasets"))

  all_metrics <- lapply(dataset_dirs, function(d) {
    tryCatch(classify_one_dataset(d), error = function(e) {
      warning("Failed dataset ", basename(d), ": ", e$message)
      NULL
    })
  })

  metrics_df <- dplyr::bind_rows(all_metrics)

  saveRDS(metrics_df, file.path(out_dir, "label_transfer_metrics.rds"))
  utils::write.csv(metrics_df, file.path(out_dir, "label_transfer_metrics.csv"), row.names = FALSE)

  writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))

  message_time("Label transfer benchmark completed")
}

if (!interactive()) {
  main()
}
