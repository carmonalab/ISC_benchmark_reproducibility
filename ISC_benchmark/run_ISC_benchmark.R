#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(yaml)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(Matrix)
  library(BiocParallel)
  library(scTypeEval)
})

source("utils/project_paths.R")
source("utils/cli_utils.R")
source("utils/benchmark_utils.R")

stop_if_not_root <- function() {
  if (!file.exists("config/benchmark_parameters.yaml") || !dir.exists("data/processed/isc")) {
    stop("Run this script from the repository root (missing config/ or data/processed/isc).")
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

load_scTypeEval_benchmarking_helpers <- function() {
  assays_utils <- system.file("benchmarking", "assays_utils.R", package = "scTypeEval")
  metrics_benchmarking <- system.file("benchmarking", "Metrics_benchmarking.R", package = "scTypeEval")
  if (assays_utils == "" || metrics_benchmarking == "") {
    stop("Could not locate scTypeEval benchmarking scripts via system.file(). ",
         "Ensure the installed scTypeEval contains inst/benchmarking/.")
  }
  source(assays_utils)
  source(metrics_benchmarking)
}

run_one_dataset <- function(dataset_path,
                            meta_path,
                            params,
                            timestamp,
                            results_root = "results/isc") {
  dataset_name <- tools::file_path_sans_ext(basename(dataset_path))

  meta <- NULL
  if (file.exists(meta_path)) {
    meta <- yaml::read_yaml(meta_path)
  }

  if (!is.null(meta$isc) && !isTRUE(meta$isc)) {
    message_time(paste("Skipping (ISC-benchmarking == no):", dataset_name))
    return(invisible(NULL))
  }

  idents <- meta$idents %||% character(0)

  obj <- scTypeEval::load_singleCell_object(dataset_path)
  counts <- obj$counts
  md <- obj$metadata

  if (!"sample" %in% colnames(md)) {
    stop("Processed dataset metadata must include a 'sample' column: ", dataset_path)
  }

  # If sidecar yaml missing, default to all ident-like columns from specs are unknown.
  if (length(idents) == 0) {
    warning("Missing sidecar yaml; defaulting idents to any character/factor columns (excluding sample)")
    candidates <- setdiff(colnames(md), c("sample", "batch", "condition", "dataset_reference", "dataset_id", "split"))
    idents <- candidates
  }

  # Keep only idents that exist.
  idents <- idents[idents %in% colnames(md)]
  if (length(idents) == 0) {
    stop("No valid idents found for dataset: ", dataset_path)
  }

  ensure_dir(results_root)

  # --- TASKS ---
  tasks_cfg <- params$isc_tasks
  tasks_out_base <- file.path(results_root, "tasks", "output")
  ensure_dir(tasks_out_base)

  task_dir <- file.path(tasks_out_base, paste0(dataset_name, "_", timestamp))
  ensure_dir(task_dir)
  yaml::write_yaml(tasks_cfg, file.path(task_dir, "params.yaml"))

  message_time(paste("ISC tasks for", dataset_name))

  for (ident in idents) {
    message_time(paste("Running tasks for ident:", ident))

    sub_meta <- md %>% filter(!is.na(.data[[ident]]))
    sub_counts <- counts[, rownames(sub_meta), drop = FALSE]

    ident_dir <- file.path(task_dir, ident)
    ensure_dir(ident_dir)

    preparams <- list(
      count_matrix = sub_counts,
      metadata = sub_meta,
      ident = ident,
      dir = ident_dir
    )

    # The wrapper functions are defined by scTypeEval's inst/benchmarking scripts.
    if (isTRUE(tasks_cfg$run$SplitCelltype)) {
      params_i <- c(preparams, tasks_cfg$common, tasks_cfg$wr.SplitCelltype)
      wr <- do.call(wr.splitCellType, params_i)
      r <- wr.assayPlot(wr, type = "Constant", return.df = TRUE) %>%
        mutate(ident = ident, dataset = dataset_name, task = "SplitCelltype")
      saveRDS(r, file.path(ident_dir, "SplitCelltype_summary.rds"))
    }

    if (isTRUE(tasks_cfg$run$missclassify)) {
      params_i <- c(preparams, tasks_cfg$common, tasks_cfg$wr.missclasify)
      wr <- do.call(wr.missclasify, params_i)
      r <- wr.assayPlot(wr, type = "monotonic", return.df = TRUE) %>%
        mutate(ident = ident, dataset = dataset_name, task = "Missclassification")
      saveRDS(r, file.path(ident_dir, "Missclassification_summary.rds"))
    }

    if (isTRUE(tasks_cfg$run$Nsamples)) {
      params_i <- c(preparams, tasks_cfg$common, tasks_cfg$wr.NSamples)
      wr <- do.call(wr.NSamples, params_i)
      r <- wr.assayPlot(wr, type = "Constant", return.df = TRUE) %>%
        mutate(ident = ident, dataset = dataset_name, task = "NSamples")
      saveRDS(r, file.path(ident_dir, "Nsamples_summary.rds"))
    }

    if (isTRUE(tasks_cfg$run$NCell)) {
      params_i <- c(preparams, tasks_cfg$common, tasks_cfg$wr.NCell)
      wr <- do.call(wr.NCell, params_i)
      r <- wr.assayPlot(wr, type = "Constant", return.df = TRUE) %>%
        mutate(ident = ident, dataset = dataset_name, task = "NCell")
      saveRDS(r, file.path(ident_dir, "NCell_summary.rds"))
    }

    if (isTRUE(tasks_cfg$run$mergeCT)) {
      params_i <- c(preparams, tasks_cfg$common, tasks_cfg$wr.mergeCT)
      wr <- do.call(wr.mergeCT, params_i)
      r <- wr.assayPlot(wr, type = "Constant", group = "ident", return.df = TRUE) %>%
        mutate(ident = ident, dataset = dataset_name, task = "mergeCT")
      saveRDS(r, file.path(ident_dir, "mergeCT_summary.rds"))
    }

    if (isTRUE(tasks_cfg$run$Nct)) {
      params_i <- c(preparams, tasks_cfg$common, tasks_cfg$wr.Nct)
      wr <- do.call(wr.Nct, params_i)
      r <- wr.assayPlot(wr, type = "Constant", return.df = TRUE) %>%
        mutate(ident = ident, dataset = dataset_name, task = "Nct")
      saveRDS(r, file.path(ident_dir, "Nct_summary.rds"))
    }
  }

  writeLines(capture.output(sessionInfo()), file.path(task_dir, "sessionInfo.txt"))

  # --- RESOURCES ---
  res_cfg <- params$isc_resources
  res_out_base <- file.path(results_root, "resources", "output")
  ensure_dir(res_out_base)

  res_dir <- file.path(res_out_base, paste0(dataset_name, "_", timestamp))
  ensure_dir(res_dir)
  yaml::write_yaml(res_cfg, file.path(res_dir, "params.yaml"))

  message_time(paste("Resource profiling for", dataset_name))

  combis <- expand.grid(res_cfg$dissimilarity.method,
                        res_cfg$IntVal.metric,
                        stringsAsFactors = FALSE)

  for (ident in idents) {
    message_time(paste("Resources for ident:", ident))

    sub_meta <- md %>% filter(!is.na(.data[[ident]]))
    sub_counts <- counts[, rownames(sub_meta), drop = FALSE]

    if (length(unique(sub_meta[[ident]])) < 2) next

    # Blacklist
    if (is.null(res_cfg$black.list)) {
      data("default_black_list", package = "scTypeEval")
      bl <- unlist(list(black.list$TCR, black.list$Immunoglobulins, black.list$Ygenes))
    } else {
      bl <- read.table(res_cfg$black.list)$V1
    }

    sc <- create.scTypeEval(matrix = sub_counts, metadata = sub_meta, black.list = bl)
    sc <- Run.ProcessingData(sc,
                             ident = ident,
                             sample = res_cfg$sample,
                             normalization.method = res_cfg$normalization.method,
                             min.samples = res_cfg$min.samples,
                             min.cells = res_cfg$min.cells,
                             verbose = res_cfg$verbose)

    if (is.null(res_cfg$gene.list)) {
      sc <- Run.HVG(sc, ncores = params$ncores %||% 1)
      gene.list <- names(sc@gene.lists[1])
    } else {
      gene.list <- res_cfg$gene.list
      sc <- add.GeneList(sc, gene.list = gene.list)
    }

    sc <- Run.PCA(sc)

    nfeatures <- length(sc@gene.lists[[1]])
    ncells <- ncol(sub_counts)
    nsamples <- length(unique(sub_meta[[res_cfg$sample]]))

    param <- BiocParallel::MulticoreParam(workers = params$ncores %||% 1,
                                          progressbar = FALSE)

    ili <- bplapply(seq_len(nrow(combis)), BPPARAM = param, function(i) {
      dm <- combis[i, 1]
      im <- combis[i, 2]

      rs <- benchmark_wrapper(expr = {
        sc.tmp <- Run.Dissimilarity(scTypeEval = sc,
                                    method = dm,
                                    reduction = res_cfg$reduction,
                                    gene.list = NULL,
                                    black.list = NULL,
                                    BestHit.classifier = "SingleR",
                                    ncores = 1,
                                    bparam = NULL,
                                    progressbar = FALSE,
                                    verbose = res_cfg$verbose)

        sc.tmp <- get.Consistency(scTypeEval = sc.tmp,
                                  dissimilarity.slot = dm,
                                  Consistency.metric = im,
                                  KNNGraph_k = res_cfg$KNNGraph_k,
                                  hclust.method = res_cfg$hclust.method,
                                  normalize = FALSE,
                                  verbose = res_cfg$verbose)
      })

      data.frame(duration = rs$duration,
                 memory_usage_MB = rs$memory / 1024^2,
                 cpu_usage = rs$cpu_usage,
                 method = im,
                 dissimilarity.method = dm)
    })

    dfi <- dplyr::bind_rows(ili) %>%
      mutate(dataset = dataset_name,
             ident = ident,
             nfeatures = nfeatures,
             ncells = ncells,
             nsamples = nsamples,
             gene.list = gene.list)

    saveRDS(dfi, file.path(res_dir, paste0(dataset_name, "_", ident, ".rds")))
  }

  writeLines(capture.output(sessionInfo()), file.path(res_dir, "sessionInfo.txt"))

  invisible(TRUE)
}

main <- function() {
  stop_if_not_root()

  load_scTypeEval_benchmarking_helpers()

  params <- yaml::read_yaml("config/benchmark_parameters.yaml")
  params <- convert_lists(params)

  timestamp <- format(Sys.time(), "%H%M%S_%d%m%Y")

  zenodo_doi <- read_zenodo_doi("config/datasets.yaml")
  if (!is.na(zenodo_doi) && nzchar(zenodo_doi)) {
    message("Zenodo DOI: ", zenodo_doi)
    ensure_dir("results/isc")
    yaml::write_yaml(
      list(
        zenodo_doi = zenodo_doi,
        run_timestamp = timestamp,
        script = "ISC_benchmark/run_ISC_benchmark.R"
      ),
      file.path("results/isc", "run_metadata.yaml")
    )
  }

  datasets <- list.files("data/processed/isc", pattern = "[.]rds$", full.names = TRUE)
  if (length(datasets) == 0) {
    stop("No processed datasets found under data/processed/isc/. Run data_processing/process_datasets.R first.")
  }

  for (ds in datasets) {
    meta_path <- paste0(ds, ".yaml")
    run_one_dataset(
      dataset_path = ds,
      meta_path = meta_path,
      params = params,
      timestamp = timestamp
    )
  }

  message_time("ISC benchmark completed")
}

if (!interactive()) {
  main()
}
