suppressPackageStartupMessages({
  library(yaml)
  library(dplyr)
  library(tidyr)
  library(Matrix)
  library(scTypeEval)
})

.resource_prepared_cache <- new.env(parent = emptyenv())

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

resource_get_prepared_input <- function(prepared_path) {
  cache_key <- normalizePath(prepared_path, mustWork = TRUE)

  if (exists(cache_key, envir = .resource_prepared_cache, inherits = FALSE)) {
    return(get(cache_key, envir = .resource_prepared_cache, inherits = FALSE))
  }

  prepared <- readRDS(cache_key)
  assign(cache_key, prepared, envir = .resource_prepared_cache)
  prepared
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

resource_time_command <- function() {
  time_cmd <- Sys.which("time")
  if (nzchar(time_cmd)) {
    return(time_cmd)
  }

  if (file.exists("/usr/bin/time")) {
    return("/usr/bin/time")
  }

  NA_character_
}

resource_parse_time_output <- function(output_lines) {
  output_lines <- output_lines[nzchar(output_lines)]
  if (length(output_lines) == 0) {
    stop("No timing output captured from /usr/bin/time")
  }

  # Use the last line that matches GNU time format: "<elapsed_seconds> <max_rss_kb>".
  candidate_idx <- grep("^[[:space:]]*[0-9]+(\\.[0-9]+)?[[:space:]]+[0-9]+[[:space:]]*$", output_lines)
  if (length(candidate_idx) == 0) {
    stop(
      "Could not find timing line in /usr/bin/time output. Output was:\n",
      paste(output_lines, collapse = "\n")
    )
  }

  timing_line <- trimws(output_lines[[tail(candidate_idx, 1)]])
  timing_values <- strsplit(timing_line, "[[:space:]]+", perl = TRUE)[[1]]
  if (length(timing_values) < 2) {
    stop("Unexpected timing output from /usr/bin/time: ", timing_line)
  }

  elapsed_seconds <- as.numeric(timing_values[[1]])
  peak_memory_kb <- as.numeric(timing_values[[2]])
  if (!is.finite(elapsed_seconds) || !is.finite(peak_memory_kb)) {
    stop("Non-numeric timing values from /usr/bin/time: ", timing_line)
  }

  list(
    elapsed_seconds = elapsed_seconds,
    peak_memory_kb = peak_memory_kb
  )
}

resource_write_benchmark_script <- function(script_path) {
  script_lines <- c(
    "project_root <- Sys.getenv('PROJECT_ROOT', unset = normalizePath(file.path(getwd(), '..')))",
    "activate <- file.path(project_root, 'renv', 'activate.R')",
    "if (file.exists(activate)) source(activate)",
    "if (requireNamespace('renv', quietly = TRUE)) {",
    "  renv::load(project = project_root)",
    "}",
    "",
    "suppressPackageStartupMessages({",
    "  library(scTypeEval)",
    "})",
    "",
    "args <- commandArgs(trailingOnly = TRUE)",
    "prepared_path <- args[[1]]",
    "result_path <- args[[2]]",
    "dissimilarity_method <- args[[3]]",
    "consistency_metric <- args[[4]]",
    "benchmark_ncores <- as.integer(args[[5]])",
    "reduction <- tolower(args[[6]]) == 'true'",
    "reciprocal_classifier <- args[[7]]",
    "knn_graph_k <- as.integer(args[[8]])",
    "hclust_method <- args[[9]]",
    "verbose_opt <- tolower(args[[10]]) == 'true'",
    "",
    "prepared <- readRDS(prepared_path)",
    "sc_tmp <- prepared$sc",
    "timing <- system.time({",
    "  sc_tmp <- scTypeEval::run_dissimilarity(",
    "    scTypeEval = sc_tmp,",
    "    method = dissimilarity_method,",
    "    reduction = reduction,",
    "    reciprocal_classifier = reciprocal_classifier,",
    "    ncores = benchmark_ncores,",
    "    verbose = verbose_opt",
    "  )",
    "",
    "  invisible(scTypeEval::get_consistency(",
    "    scTypeEval = sc_tmp,",
    "    dissimilarity_slot = dissimilarity_method,",
    "    consistency_metric = consistency_metric,",
    "    knn_graph_k = knn_graph_k,",
    "    hclust_method = hclust_method,",
    "    normalize = FALSE,",
    "    verbose = verbose_opt",
    "  ))",
    "})",
    "",
    "duration_ms <- unname(timing[['elapsed']]) * 1000",
    "cpu_usage <- if (timing[['elapsed']] > 0) {",
    "  (timing[['user.self']] + timing[['sys.self']]) / timing[['elapsed']]",
    "} else {",
    "  NA_real_",
    "}",
    "",
    "saveRDS(",
    "  list(duration_ms = duration_ms, cpu_usage = cpu_usage),",
    "  result_path",
    ")"
  )

  writeLines(script_lines, script_path)
  script_path
}

resource_run_single_benchmark <- function(prepared_path,
                                          dissimilarity_method,
                                          consistency_metric,
                                          params) {
  time_cmd <- resource_time_command()
  if (is.na(time_cmd)) {
    stop("/usr/bin/time not found; peak memory reporting is unavailable on this system")
  }

  benchmark_ncores <- 1L
  benchmark_iterations <- as.integer(params$benchmark$iterations)
  if (is.na(benchmark_iterations) || benchmark_iterations < 1L) {
    benchmark_iterations <- 1L
  }

  configured_cores <- as.integer(params$benchmark$ncores)
  if (!is.na(configured_cores) && configured_cores != 1L) {
    warning(
      "Overriding benchmark.ncores=", configured_cores,
      " to 1 for resource measurement.",
      call. = FALSE
    )
  }

  child_script <- tempfile("resource_benchmark_", fileext = ".R")
  child_result <- tempfile("resource_benchmark_result_", fileext = ".rds")
  on.exit(unlink(c(child_script, child_result)), add = TRUE)

  resource_write_benchmark_script(child_script)

  if (!nzchar(Sys.which("Rscript"))) {
    stop("Rscript not found on PATH")
  }

  duration_ms_values <- numeric(benchmark_iterations)
  cpu_usage_values <- numeric(benchmark_iterations)
  peak_memory_mb_values <- numeric(benchmark_iterations)
  elapsed_seconds_values <- numeric(benchmark_iterations)

  for (i in seq_len(benchmark_iterations)) {
    command_args <- c(
      shQuote(time_cmd),
      "-f",
      shQuote("%e\t%M"),
      shQuote(Sys.which("Rscript")),
      "--vanilla",
      shQuote(child_script),
      shQuote(prepared_path),
      shQuote(child_result),
      shQuote(dissimilarity_method),
      shQuote(consistency_metric),
      shQuote(as.character(benchmark_ncores)),
      shQuote(as.character(isTRUE(params$common$reduction))),
      shQuote(params$common$reciprocal_classifier),
      shQuote(as.character(params$common$knn_graph_k)),
      shQuote(params$common$hclust_method),
      shQuote(as.character(isTRUE(params$common$verbose)))
    )

    shell_command <- paste(command_args, collapse = " ")
    shell_command <- paste(shell_command, "2>&1")

    timing_output <- system(
      shell_command,
      intern = TRUE,
      ignore.stderr = FALSE
    )

    exit_status <- attr(timing_output, "status")
    if (!is.null(exit_status) && !identical(exit_status, 0L)) {
      stop(
        "Resource benchmark failed for ", dissimilarity_method, " / ", consistency_metric,
        " at replicate ", i, " of ", benchmark_iterations,
        " with exit status ", exit_status, ". Output:\n",
        paste(timing_output, collapse = "\n")
      )
    }

    timing <- resource_parse_time_output(timing_output)
    benchmark_summary <- readRDS(child_result)

    duration_ms_values[[i]] <- as.numeric(benchmark_summary$duration_ms)
    cpu_usage_values[[i]] <- as.numeric(benchmark_summary$cpu_usage)
    peak_memory_mb_values[[i]] <- as.numeric(timing$peak_memory_kb) / 1024
    elapsed_seconds_values[[i]] <- as.numeric(timing$elapsed_seconds)
  }

  list(
    duration_ms = stats::median(duration_ms_values, na.rm = TRUE),
    cpu_usage = stats::median(cpu_usage_values, na.rm = TRUE),
    peak_memory_MB = stats::median(peak_memory_mb_values, na.rm = TRUE),
    benchmark_ncores = benchmark_ncores,
    benchmark_iterations = benchmark_iterations,
    benchmark_elapsed_seconds = stats::median(elapsed_seconds_values, na.rm = TRUE)
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
  prepared <- resource_get_prepared_input(prepared_path)

  resource_message_time(
    "Benchmarking ", prepared$dataset_id, " / ", prepared$ident,
    " / ", consistency_metric, " / ", dissimilarity_method
  )

  benchmark_summary <- resource_run_single_benchmark(
    prepared_path = prepared_path,
    dissimilarity_method = dissimilarity_method,
    consistency_metric = consistency_metric,
    params = params
  )

  output <- data.frame(
    duration_ms = benchmark_summary$duration_ms,
    peak_memory_MB = benchmark_summary$peak_memory_MB,
    duration = benchmark_summary$duration_ms,
    memory_usage_MB = benchmark_summary$peak_memory_MB,
    cpu_usage = benchmark_summary$cpu_usage,
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
    benchmark_ncores = benchmark_summary$benchmark_ncores,
    benchmark_iterations = benchmark_summary$benchmark_iterations,
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

# Check if a dataset/ident combination has been processed
is_dataset_ident_completed <- function(params, dataset_id, ident) {
  output_root <- resource_output_dir_from_config(params)
  dataset_dir <- file.path(output_root, sanitize_for_path(dataset_id))
  ident_dir <- file.path(dataset_dir, sanitize_for_path(ident))

  if (!dir.exists(ident_dir)) {
    return(FALSE)
  }

  expected_n <- length(params$common$dissimilarity_method) * length(params$common$consistency_metric)
  result_files <- list.files(ident_dir, pattern = ".*\\.rds$", full.names = TRUE)

  if (length(result_files) < expected_n) {
    return(FALSE)
  }

  # Ensure previously produced outputs contain usable memory metrics.
  is_valid_result <- function(path) {
    x <- tryCatch(readRDS(path), error = function(e) NULL)
    if (is.null(x)) {
      return(FALSE)
    }

    peak <- if ("peak_memory_MB" %in% names(x)) x$peak_memory_MB else x$memory_usage_MB
    dur <- if ("duration_ms" %in% names(x)) x$duration_ms else x$duration
    expected_iterations <- as.integer(params$benchmark$iterations)
    if (is.na(expected_iterations) || expected_iterations < 1L) {
      expected_iterations <- 1L
    }
    result_iterations <- if ("benchmark_iterations" %in% names(x)) {
      as.integer(x$benchmark_iterations)
    } else {
      NA_integer_
    }

    is.finite(as.numeric(peak)) &&
      is.finite(as.numeric(dur)) &&
      !is.na(result_iterations) &&
      identical(result_iterations, expected_iterations)
  }

  all(vapply(result_files, is_valid_result, logical(1)))
}

# Filter to only incomplete dataset/ident combinations
filter_incomplete_ident_grid <- function(ident_grid, params) {
  if (is.null(ident_grid) || nrow(ident_grid) == 0) {
    return(ident_grid)
  }

  incomplete <- mapply(
    function(dataset_id, ident) {
      !is_dataset_ident_completed(params, dataset_id, ident)
    },
    ident_grid$dataset_id,
    ident_grid$ident,
    SIMPLIFY = TRUE
  )

  ident_grid[incomplete, , drop = FALSE]
}