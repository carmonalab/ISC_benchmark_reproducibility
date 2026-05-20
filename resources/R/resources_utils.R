suppressPackageStartupMessages({
  library(yaml)
  library(dplyr)
  library(tidyr)
  library(Matrix)
  library(scTypeEval)
})

resource_proj_root <- function(start_dir = getwd()) {
  current <- normalizePath(start_dir, mustWork = TRUE)
  max_depth <- 10
  depth <- 0

  while (depth < max_depth) {
    if (length(list.files(current, pattern = "\\.Rproj$")) > 0 ||
        file.exists(file.path(current, "renv.lock"))) {
      return(current)
    }

    parent <- dirname(current)
    if (identical(parent, current)) {
      break
    }

    current <- parent
    depth <- depth + 1
  }

  stop("Could not determine project root for resources pipeline")
}

resource_proj_path <- function(...) {
  file.path(resource_proj_root(), ...)
}

resource_ensure_dir <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
  invisible(path)
}

resource_convert_lists <- function(x) {
  if (is.list(x)) {
    x <- lapply(x, resource_convert_lists)
    if (is.null(names(x)) && all(vapply(x, is.atomic, logical(1)))) {
      return(unlist(x, use.names = FALSE))
    }
  }
  x
}

resource_message_time <- function(...) {
  msg <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", ...)
  message(msg)
}

sanitize_for_path <- function(x) {
  gsub("[^A-Za-z0-9._-]+", "_", x)
}

resource_output_dir <- function() {
  resource_ensure_dir(resource_proj_path("resources", "output"))
}

resource_cache_dir <- function() {
  resource_ensure_dir(resource_proj_path("resources", "cache"))
}

resource_prepared_dir <- function() {
  resource_ensure_dir(file.path(resource_cache_dir(), "prepared"))
}

resource_output_dir_from_config <- function(params) {
  resource_ensure_dir(params$output_root)
}

resource_cache_dir_from_config <- function(params) {
  resource_ensure_dir(params$cache_root)
}

resource_prepared_dir_from_config <- function(params) {
  resource_ensure_dir(file.path(resource_cache_dir_from_config(params), "prepared"))
}

resolve_root_path <- function(path_value, mustWork = TRUE) {
  if (startsWith(path_value, "/")) {
    return(normalizePath(path_value, mustWork = mustWork))
  }
  normalizePath(resource_proj_path(path_value), mustWork = mustWork)
}

load_resource_config <- function() {
  config <- yaml::read_yaml("config/resource_parameters.yaml")
  config <- resource_convert_lists(config)

  config$processed_data_dir <- resolve_root_path(config$paths$data_processed, mustWork = TRUE)
  config$dataset_idents_file <- resolve_root_path(config$paths$dataset_idents_file, mustWork = TRUE)
  config$output_root <- resolve_root_path(config$paths$output_root, mustWork = FALSE)
  config$cache_root <- resolve_root_path(config$paths$cache_root, mustWork = FALSE)

  resource_ensure_dir(config$output_root)
  resource_ensure_dir(config$cache_root)
  resource_ensure_dir(file.path(config$cache_root, "prepared"))

  if (!is.null(config$common$black_list)) {
    config$common$black_list <- resolve_root_path(config$common$black_list, mustWork = TRUE)
  }
  if (!is.null(config$common$gene_list)) {
    config$common$gene_list <- resolve_root_path(config$common$gene_list, mustWork = TRUE)
  }

  config
}

get_requested_resource_dataset_ids <- function() {
  requested <- c(
    Sys.getenv("RESOURCE_DATASET_ID", unset = ""),
    Sys.getenv("RESOURCE_TEST_DATASET", unset = ""),
    Sys.getenv("RESOURCE_DATASET_IDS", unset = "")
  )

  requested <- requested[nzchar(requested)]
  if (length(requested) == 0) {
    return(NULL)
  }

  unique(trimws(unlist(strsplit(paste(requested, collapse = ","), ",", fixed = TRUE))))
}

normalize_dataset_family <- function(dataset_id) {
  sub("_.*$", "", dataset_id)
}

get_resource_dataset_idents <- function(config, selected_dataset_ids = NULL) {
  ident_mapping <- yaml::read_yaml(config$dataset_idents_file)$idents
  files <- list.files(config$processed_data_dir, pattern = "\\.rds$", full.names = TRUE)

  if (length(files) == 0) {
    stop("No processed datasets found in: ", config$processed_data_dir)
  }

  dataset_info <- lapply(files, function(dataset_file) {
    dataset_id <- tools::file_path_sans_ext(basename(dataset_file))
    dataset_family <- normalize_dataset_family(dataset_id)
    ident_cols <- ident_mapping[[dataset_family]]

    if (is.null(ident_cols) || length(ident_cols) == 0) {
      warning("No ident columns configured for dataset family: ", dataset_family)
      return(NULL)
    }

    data.frame(
      dataset_id = dataset_id,
      dataset_file = dataset_file,
      dataset_family = dataset_family,
      ident_cols = paste(ident_cols, collapse = ","),
      stringsAsFactors = FALSE
    )
  })

  dataset_info <- Filter(Negate(is.null), dataset_info)
  if (length(dataset_info) == 0) {
    stop("No valid processed datasets matched dataset_idents configuration")
  }

  dataset_info <- dplyr::bind_rows(dataset_info)

  if (!is.null(selected_dataset_ids)) {
    selected_dataset_ids <- unique(trimws(selected_dataset_ids))
    selected_dataset_ids <- selected_dataset_ids[nzchar(selected_dataset_ids)]

    missing_ids <- setdiff(selected_dataset_ids, dataset_info$dataset_id)
    if (length(missing_ids) > 0) {
      stop(
        "Requested dataset(s) not available in processed inputs: ",
        paste(missing_ids, collapse = ", ")
      )
    }

    dataset_info <- dataset_info[match(selected_dataset_ids, dataset_info$dataset_id), , drop = FALSE]
  }

  dataset_info
}

build_resource_ident_grid <- function(dataset_info) {
  if (is.null(dataset_info) || nrow(dataset_info) == 0) {
    stop("No dataset information available to build resource grid")
  }

  rows <- lapply(seq_len(nrow(dataset_info)), function(i) {
    ident_cols <- trimws(unlist(strsplit(dataset_info$ident_cols[[i]], ",", fixed = TRUE)))
    ident_cols <- ident_cols[nzchar(ident_cols)]

    data.frame(
      dataset_id = dataset_info$dataset_id[[i]],
      dataset_file = dataset_info$dataset_file[[i]],
      ident = ident_cols,
      stringsAsFactors = FALSE
    )
  })

  dplyr::bind_rows(rows)
}

build_resource_metric_grid <- function(config) {
  tidyr::expand_grid(
    dissimilarity_method = config$common$dissimilarity_method,
    consistency_metric = config$common$consistency_metric
  )
}

read_optional_gene_list <- function(path) {
  if (is.null(path)) {
    return(NULL)
  }

  gene_values <- read.table(path, stringsAsFactors = FALSE)[[1]]
  gene_values <- unique(gene_values[nzchar(gene_values)])
  list(custom = gene_values)
}

benchmark_wrapper <- function(expr, iterations = 3, envir = parent.frame()) {
  durations <- numeric(iterations)
  cpu_percent <- numeric(iterations)
  expr <- base::substitute(expr)

  for (i in seq_len(iterations)) {
    gc(verbose = FALSE)
    timing <- system.time(eval(expr, envir = envir))
    durations[[i]] <- unname(timing[["elapsed"]])
    cpu_percent[[i]] <- if (timing[["elapsed"]] > 0) {
      (timing[["user.self"]] + timing[["sys.self"]]) / timing[["elapsed"]]
    } else {
      NA_real_
    }
  }

  gc(verbose = FALSE)
  memory_profile <- tempfile(fileext = ".out")
  on.exit(unlink(memory_profile), add = TRUE)

  Rprofmem(memory_profile)
  on.exit(Rprofmem(NULL), add = TRUE)
  eval(expr, envir = envir)
  Rprofmem(NULL)

  profile_lines <- readLines(memory_profile, warn = FALSE)
  memory_values <- suppressWarnings(as.numeric(sub(" .*", "", profile_lines)))

  data.frame(
    duration = stats::median(durations, na.rm = TRUE),
    cpu_usage = stats::median(cpu_percent, na.rm = TRUE),
    memory = sum(memory_values[is.finite(memory_values)], na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

prepare_resource_input <- function(dataset_id, dataset_file, ident, params) {
  resource_message_time("Preparing resource benchmark input for ", dataset_id, " / ", ident)

  loaded <- scTypeEval::load_single_cell_object(dataset_file)
  metadata <- as.data.frame(loaded$metadata)
  sample_col <- params$common$sample

  if (!ident %in% colnames(metadata)) {
    stop("Cell type column not found in metadata: ", ident)
  }
  if (!sample_col %in% colnames(metadata)) {
    stop("Sample column not found in metadata: ", sample_col)
  }

  valid_cells <- !is.na(metadata[[ident]]) & !is.na(metadata[[sample_col]])
  metadata <- metadata[valid_cells, , drop = FALSE]
  count_matrix <- loaded$counts[, rownames(metadata), drop = FALSE]

  black_list <- NULL
  if (!is.null(params$common$black_list)) {
    black_list <- read.table(params$common$black_list, stringsAsFactors = FALSE)[[1]]
  }

  sc <- scTypeEval::create_scTypeEval(
    matrix = count_matrix,
    metadata = metadata,
    active_ident = ident,
    black_list = black_list
  )

  sc <- scTypeEval::run_processing_data(
    scTypeEval = sc,
    ident = ident,
    sample = sample_col,
    normalization_method = params$common$normalization_method,
    min_samples = params$common$min_samples,
    min_cells = params$common$min_cells,
    verbose = isTRUE(params$common$verbose)
  )

  gene_list <- read_optional_gene_list(params$common$gene_list)
  if (is.null(gene_list)) {
    sc <- scTypeEval::run_hvg(
      scTypeEval = sc,
      ncores = params$n_cores,
      verbose = isTRUE(params$common$verbose)
    )
  } else {
    sc <- scTypeEval::add_gene_list(scTypeEval = sc, gene_list = gene_list)
  }

  if (isTRUE(params$common$reduction)) {
    sc <- scTypeEval::run_pca(
      scTypeEval = sc,
      ndim = params$common$ndim,
      verbose = isTRUE(params$common$verbose)
    )
  }

  gene_list_name <- names(sc@gene_lists)[[1]]
  feature_set <- intersect(sc@gene_lists[[gene_list_name]], rownames(count_matrix))
  nfeatures <- length(feature_set)

  feature_matrix <- count_matrix[feature_set, , drop = FALSE]
  total_entries <- prod(dim(feature_matrix))
  sparsity <- if (total_entries == 0) {
    NA_real_
  } else {
    (total_entries - Matrix::nnzero(feature_matrix)) / total_entries
  }

  prepared_path <- file.path(
    resource_prepared_dir_from_config(params),
    paste0(sanitize_for_path(dataset_id), "__", sanitize_for_path(ident), ".rds")
  )

  saveRDS(
    list(
      dataset_id = dataset_id,
      dataset_file = dataset_file,
      ident = ident,
      sample_col = sample_col,
      sc = sc,
      nfeatures = nfeatures,
      ncells = ncol(count_matrix),
      nsamples = dplyr::n_distinct(metadata[[sample_col]]),
      sparsity = sparsity,
      gene_list_name = gene_list_name
    ),
    prepared_path
  )

  prepared_path
}

build_resource_output_path <- function(params,
                                       dataset_id,
                                       ident,
                                       consistency_metric,
                                       dissimilarity_method) {
  dataset_dir <- resource_ensure_dir(file.path(resource_output_dir_from_config(params), sanitize_for_path(dataset_id)))
  ident_dir <- resource_ensure_dir(file.path(dataset_dir, sanitize_for_path(ident)))

  file.path(
    ident_dir,
    paste0(
      sanitize_for_path(consistency_metric),
      "__",
      sanitize_for_path(dissimilarity_method),
      ".rds"
    )
  )
}

benchmark_resource_pair <- function(prepared_path, dissimilarity_method, consistency_metric, params) {
  prepared <- readRDS(prepared_path)
  sc_template <- prepared$sc
  iterations <- as.integer(params$benchmark$iterations)
  benchmark_ncores <- as.integer(params$benchmark$ncores)
  verbose_opt <- isTRUE(params$common$verbose)

  resource_message_time(
    "Benchmarking ", prepared$dataset_id, " / ", prepared$ident,
    " / ", consistency_metric, " / ", dissimilarity_method
  )

  stats <- benchmark_wrapper(
    expr = {
      sc_tmp <- sc_template
      sc_tmp <- scTypeEval::run_dissimilarity(
        scTypeEval = sc_tmp,
        method = dissimilarity_method,
        reduction = isTRUE(params$common$reduction),
        reciprocal_classifier = params$common$reciprocal_classifier,
        ncores = benchmark_ncores,
        verbose = verbose_opt
      )

      invisible(scTypeEval::get_consistency(
        scTypeEval = sc_tmp,
        dissimilarity_slot = dissimilarity_method,
        consistency_metric = consistency_metric,
        knn_graph_k = params$common$knn_graph_k,
        hclust_method = params$common$hclust_method,
        normalize = FALSE,
        verbose = verbose_opt
      ))
    },
    iterations = iterations
  )

  output <- data.frame(
    duration = stats$duration,
    memory_usage_MB = stats$memory / 1024^2,
    cpu_usage = stats$cpu_usage,
    method = consistency_metric,
    consistency_metric = consistency_metric,
    consistency.metric = consistency_metric,
    dissimilarity_method = dissimilarity_method,
    dissimilarity.method = dissimilarity_method,
    dataset = prepared$dataset_id,
    dataset_id = prepared$dataset_id,
    ident = prepared$ident,
    nfeatures = prepared$nfeatures,
    ncells = prepared$ncells,
    nsamples = prepared$nsamples,
    sparsity = prepared$sparsity,
    gene.list = prepared$gene_list_name,
    benchmark_ncores = benchmark_ncores,
    benchmark_iterations = iterations,
    stringsAsFactors = FALSE
  )

  output_path <- build_resource_output_path(
    params = params,
    dataset_id = prepared$dataset_id,
    ident = prepared$ident,
    consistency_metric = consistency_metric,
    dissimilarity_method = dissimilarity_method
  )

  saveRDS(output, output_path)
  output_path
}