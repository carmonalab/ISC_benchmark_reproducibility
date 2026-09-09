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

# Unified grid of every internal (scTypeEval) dissimilarity/int_val_metric pair
# plus every enabled external tool (SCCAF, anticor_features, sc-SHC), so both
# R and Python tools flow through the same benchmarking/aggregation path.
build_resource_tool_grid <- function(config) {
  internal_grid <- tidyr::expand_grid(
    dissimilarity_method = config$common$dissimilarity_method,
    int_val_metric = config$common$int_val_metric
  )
  internal_grid$tool_type <- "internal"
  internal_grid$tool_name <- paste(internal_grid$dissimilarity_method, internal_grid$int_val_metric, sep = "::")
  internal_grid$language <- "R"

  ext_cfg <- config$external_methods
  external_rows <- list()

  if (!is.null(ext_cfg) && isTRUE(ext_cfg$enabled)) {
    if (!is.null(ext_cfg$sccaf) && isTRUE(ext_cfg$sccaf$enabled)) {
      external_rows[["sccaf"]] <- data.frame(
        dissimilarity_method = NA_character_, int_val_metric = NA_character_,
        tool_type = "external_py", tool_name = "sccaf", language = "python",
        stringsAsFactors = FALSE
      )
    }
    if (!is.null(ext_cfg$anticor_features) && isTRUE(ext_cfg$anticor_features$enabled)) {
      external_rows[["anticor_features"]] <- data.frame(
        dissimilarity_method = NA_character_, int_val_metric = NA_character_,
        tool_type = "external_py", tool_name = "anticor_features", language = "python",
        stringsAsFactors = FALSE
      )
    }
    if (!is.null(ext_cfg$scshc) && isTRUE(ext_cfg$scshc$enabled)) {
      external_rows[["scshc"]] <- data.frame(
        dissimilarity_method = NA_character_, int_val_metric = NA_character_,
        tool_type = "external_r", tool_name = "scshc", language = "R",
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(external_rows) > 0) {
    dplyr::bind_rows(internal_grid, dplyr::bind_rows(external_rows))
  } else {
    internal_grid
  }
}

# Prefix that pins BLAS/OMP thread pools to `ncores` so R and Python workloads
# are measured under the same single-threaded constraint.
resource_thread_limit_prefix <- function(ncores = 1L) {
  vars <- c("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS")
  paste(paste0(vars, "=", as.integer(ncores)), collapse = " ")
}

# Lazily source the Python-pipeline helpers (pipeline_spec_*, .find_repo_root,
# .default_python_bin, .resolve_script_path, .export_sctypeeval_to_h5ad) by
# absolute path so this works regardless of the caller's working directory.
resource_load_python_pipeline_helpers <- function() {
  flag <- ".resource_python_helpers_loaded"
  if (!exists(flag, envir = .resource_prepared_cache, inherits = FALSE)) {
    source(resource_proj_path("sample_agnostic_utils", "run_python_pipelines.R"))
    assign(flag, TRUE, envir = .resource_prepared_cache)
  }
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

  # GNU time format: "<elapsed_seconds> <max_rss_kb> <user_seconds> <sys_seconds>".
  numeric_field <- "[0-9]+(\\.[0-9]+)?"
  pattern <- paste0(
    "^[[:space:]]*", numeric_field, "[[:space:]]+[0-9]+[[:space:]]+",
    numeric_field, "[[:space:]]+", numeric_field, "[[:space:]]*$"
  )
  candidate_idx <- grep(pattern, output_lines)
  if (length(candidate_idx) == 0) {
    stop(
      "Could not find timing line in /usr/bin/time output. Output was:\n",
      paste(output_lines, collapse = "\n")
    )
  }

  timing_line <- trimws(output_lines[[tail(candidate_idx, 1)]])
  timing_values <- strsplit(timing_line, "[[:space:]]+", perl = TRUE)[[1]]
  if (length(timing_values) < 4) {
    stop("Unexpected timing output from /usr/bin/time: ", timing_line)
  }

  elapsed_seconds <- as.numeric(timing_values[[1]])
  peak_memory_kb <- as.numeric(timing_values[[2]])
  user_seconds <- as.numeric(timing_values[[3]])
  sys_seconds <- as.numeric(timing_values[[4]])
  if (!all(is.finite(c(elapsed_seconds, peak_memory_kb, user_seconds, sys_seconds)))) {
    stop("Non-numeric timing values from /usr/bin/time: ", timing_line)
  }

  list(
    elapsed_seconds = elapsed_seconds,
    peak_memory_kb = peak_memory_kb,
    user_seconds = user_seconds,
    sys_seconds = sys_seconds
  )
}

resource_write_internal_benchmark_script <- function(script_path) {
  script_lines <- c(
    "project_root <- Sys.getenv('PROJECT_ROOT', unset = normalizePath(file.path(getwd(), '..')))",
    "Sys.setenv(RENV_PROJECT_EXPLICIT = project_root)",
    "Sys.setenv(RENV_CONFIG_AUTOLOADER_ENABLED = 'FALSE')",
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
    "dissimilarity_method <- args[[2]]",
    "int_val_metric <- args[[3]]",
    "benchmark_ncores <- as.integer(args[[4]])",
    "reduction <- tolower(args[[5]]) == 'true'",
    "reciprocal_classifier <- args[[6]]",
    "knn_graph_k <- as.integer(args[[7]])",
    "hclust_method <- args[[8]]",
    "verbose_opt <- tolower(args[[9]]) == 'true'",
    "",
    "prepared <- readRDS(prepared_path)",
    "sc_tmp <- prepared$sc",
    "sc_tmp <- scTypeEval::run_dissimilarity(",
    "  scTypeEval = sc_tmp,",
    "  method = dissimilarity_method,",
    "  reduction = reduction,",
    "  reciprocal_classifier = reciprocal_classifier,",
    "  ncores = benchmark_ncores,",
    "  verbose = verbose_opt",
    ")",
    "",
    "invisible(scTypeEval::get_consistency(",
    "  scTypeEval = sc_tmp,",
    "  dissimilarity_slot = dissimilarity_method,",
    "  consistency_metric = int_val_metric,",
    "  knn_graph_k = knn_graph_k,",
    "  hclust_method = hclust_method,",
    "  normalize = FALSE,",
    "  verbose = verbose_opt",
    "))"
  )

  writeLines(script_lines, script_path)
  script_path
}

resource_write_scshc_benchmark_script <- function(script_path) {
  script_lines <- c(
    "project_root <- Sys.getenv('PROJECT_ROOT', unset = normalizePath(file.path(getwd(), '..')))",
    "Sys.setenv(RENV_PROJECT_EXPLICIT = project_root)",
    "Sys.setenv(RENV_CONFIG_AUTOLOADER_ENABLED = 'FALSE')",
    "activate <- file.path(project_root, 'renv', 'activate.R')",
    "if (file.exists(activate)) source(activate)",
    "if (requireNamespace('renv', quietly = TRUE)) {",
    "  renv::load(project = project_root)",
    "}",
    "",
    "suppressPackageStartupMessages({",
    "  library(scTypeEval)",
    "})",
    "source(file.path(project_root, 'sample_agnostic_utils', 'scshc.R'))",
    "",
    "args <- commandArgs(trailingOnly = TRUE)",
    "prepared_path <- args[[1]]",
    "benchmark_ncores <- as.integer(args[[2]])",
    "",
    "prepared <- readRDS(prepared_path)",
    "invisible(run_scSHC(",
    "  scTypeEval = prepared$sc,",
    "  parallel = benchmark_ncores > 1,",
    "  cores = benchmark_ncores",
    "))"
  )

  writeLines(script_lines, script_path)
  script_path
}

# Build the (executable, args) for benchmarking an external Python tool
# directly (no R wrapper process), so timings reflect only the tool itself.
resource_build_external_python_command <- function(tool_name, h5ad_path, output_csv, params) {
  resource_load_python_pipeline_helpers()
  repo_root <- .find_repo_root()
  ext_cfg <- params$external_methods

  spec <- switch(
    tool_name,
    sccaf = pipeline_spec_sccaf(
      n = as.integer(ext_cfg$sccaf$params$n %||% 100)
    ),
    anticor_features = pipeline_spec_anticor_features(
      min_cells = as.integer(ext_cfg$anticor_features$params$min_cells %||% 10),
      species = ext_cfg$anticor_features$params$species %||% "hsapiens",
      score_k = as.numeric(ext_cfg$anticor_features$params$score_k %||% 1.0)
    ),
    stop("Unsupported external python tool: ", tool_name)
  )

  script_path <- .resolve_script_path(spec$script, repo_root)
  python_bin <- if (identical(tool_name, "sccaf")) {
    .default_sccaf_python_bin(repo_root)
  } else {
    .default_python_bin(repo_root)
  }

  if (!file.exists(python_bin)) {
    stop("Python interpreter not found for ", tool_name, ": ", python_bin)
  }

  list(
    executable = python_bin,
    args = c(script_path, .input_flag_for_script(script_path), h5ad_path, "--output", output_csv, spec$args)
  )
}

# Run `executable args...` under `/usr/bin/time`, `benchmark_iterations` times,
# and return the median duration/CPU/peak-memory. Used identically for
# internal (R) and external (R or Python) tools so metrics are comparable.
resource_run_measured_command <- function(executable, args, iterations, ncores = 1L) {
  time_cmd <- resource_time_command()
  if (is.na(time_cmd)) {
    stop("/usr/bin/time not found; peak memory reporting is unavailable on this system")
  }

  thread_env_prefix <- resource_thread_limit_prefix(ncores)

  duration_ms_values <- numeric(iterations)
  cpu_usage_values <- numeric(iterations)
  peak_memory_mb_values <- numeric(iterations)

  for (i in seq_len(iterations)) {
    quoted_command <- vapply(c(executable, args), shQuote, character(1))
    shell_command <- paste(
      thread_env_prefix,
      shQuote(time_cmd), "-f", shQuote("%e\t%M\t%U\t%S"),
      paste(quoted_command, collapse = " "),
      "2>&1"
    )

    timing_output <- system(shell_command, intern = TRUE, ignore.stderr = FALSE)

    exit_status <- attr(timing_output, "status")
    if (!is.null(exit_status) && !identical(exit_status, 0L)) {
      stop(
        "Resource benchmark command failed at replicate ", i, " of ", iterations,
        " for ", basename(executable), " with exit status ", exit_status, ". Output:\n",
        paste(timing_output, collapse = "\n")
      )
    }

    timing <- resource_parse_time_output(timing_output)
    duration_ms_values[[i]] <- timing$elapsed_seconds * 1000
    cpu_usage_values[[i]] <- if (timing$elapsed_seconds > 0) {
      (timing$user_seconds + timing$sys_seconds) / timing$elapsed_seconds
    } else {
      NA_real_
    }
    peak_memory_mb_values[[i]] <- timing$peak_memory_kb / 1024
  }

  list(
    duration_ms = stats::median(duration_ms_values, na.rm = TRUE),
    cpu_usage = stats::median(cpu_usage_values, na.rm = TRUE),
    peak_memory_MB = stats::median(peak_memory_mb_values, na.rm = TRUE),
    benchmark_ncores = as.integer(ncores),
    benchmark_iterations = as.integer(iterations)
  )
}

# Dispatch to the right measured command for a tool_grid row (internal
# scTypeEval pair, external R tool, or external Python tool).
resource_run_tool_benchmark <- function(prepared, tool_row, params) {
  iterations <- as.integer(params$benchmark$iterations)
  if (is.na(iterations) || iterations < 1L) {
    iterations <- 1L
  }

  ncores <- 1L
  configured_cores <- as.integer(params$benchmark$ncores)
  if (!is.na(configured_cores) && configured_cores != 1L) {
    warning(
      "Overriding benchmark.ncores=", configured_cores,
      " to 1 for resource measurement.",
      call. = FALSE
    )
  }

  if (identical(tool_row$tool_type, "internal")) {
    if (!nzchar(Sys.which("Rscript"))) stop("Rscript not found on PATH")
    child_script <- tempfile("resource_benchmark_", fileext = ".R")
    on.exit(unlink(child_script), add = TRUE)
    resource_write_internal_benchmark_script(child_script)

    resource_run_measured_command(
      executable = Sys.which("Rscript"),
      args = c(
        "--vanilla", child_script,
        prepared$prepared_path,
        tool_row$dissimilarity_method,
        tool_row$int_val_metric,
        as.character(ncores),
        as.character(isTRUE(params$common$reduction)),
        params$common$reciprocal_classifier,
        as.character(params$common$knn_graph_k),
        params$common$hclust_method,
        as.character(isTRUE(params$common$verbose))
      ),
      iterations = iterations,
      ncores = ncores
    )
  } else if (identical(tool_row$tool_type, "external_r")) {
    if (!nzchar(Sys.which("Rscript"))) stop("Rscript not found on PATH")
    child_script <- tempfile("resource_benchmark_scshc_", fileext = ".R")
    on.exit(unlink(child_script), add = TRUE)
    resource_write_scshc_benchmark_script(child_script)

    resource_run_measured_command(
      executable = Sys.which("Rscript"),
      args = c("--vanilla", child_script, prepared$prepared_path, as.character(ncores)),
      iterations = iterations,
      ncores = ncores
    )
  } else if (identical(tool_row$tool_type, "external_py")) {
    if (is.null(prepared$h5ad_path) || !file.exists(prepared$h5ad_path)) {
      stop("h5ad export missing for external python benchmarking of ", prepared$dataset_id)
    }
    output_csv <- tempfile(paste0("resource_benchmark_", tool_row$tool_name, "_"), fileext = ".csv")
    on.exit(unlink(output_csv), add = TRUE)

    command <- resource_build_external_python_command(
      tool_name = tool_row$tool_name,
      h5ad_path = prepared$h5ad_path,
      output_csv = output_csv,
      params = params
    )

    resource_run_measured_command(
      executable = command$executable,
      args = command$args,
      iterations = iterations,
      ncores = ncores
    )
  } else {
    stop("Unsupported tool_type: ", tool_row$tool_type)
  }
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

  # Export an h5ad alongside the prepared object when any external Python
  # tool is enabled, so external_py benchmarks measure only the tool itself.
  ext_cfg <- params$external_methods
  external_py_enabled <- !is.null(ext_cfg) && isTRUE(ext_cfg$enabled) && (
    (!is.null(ext_cfg$sccaf) && isTRUE(ext_cfg$sccaf$enabled)) ||
      (!is.null(ext_cfg$anticor_features) && isTRUE(ext_cfg$anticor_features$enabled))
  )

  h5ad_path <- NULL
  if (external_py_enabled) {
    resource_load_python_pipeline_helpers()
    h5ad_path <- file.path(
      resource_prepared_dir_from_config(params),
      paste0(sanitize_for_path(dataset_id), "__", sanitize_for_path(ident), ".h5ad")
    )
    .export_sctypeeval_to_h5ad(sc, h5ad_path)
  }

  saveRDS(
    list(
      dataset_id = dataset_id,
      dataset_file = dataset_file,
      ident = ident,
      sample_col = sample_col,
      sc = sc,
      prepared_path = prepared_path,
      h5ad_path = h5ad_path,
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

build_resource_output_path <- function(params, dataset_id, ident, tool_name) {
  dataset_dir <- resource_ensure_dir(file.path(resource_output_dir_from_config(params), sanitize_for_path(dataset_id)))
  ident_dir <- resource_ensure_dir(file.path(dataset_dir, sanitize_for_path(ident)))

  file.path(ident_dir, paste0(sanitize_for_path(tool_name), ".rds"))
}

# Benchmark a single row of build_resource_tool_grid() (internal ISC pair or
# external tool) for one prepared dataset/ident. Output schema is shared by
# internal and external rows (plus legacy aliases kept for downstream Rmds).
benchmark_resource_tool <- function(prepared_path, dissimilarity_method, int_val_metric, tool_type, tool_name, language, params) {
  prepared <- resource_get_prepared_input(prepared_path)

  resource_message_time(
    "Benchmarking ", prepared$dataset_id, " / ", prepared$ident, " / ", tool_name
  )

  tool_row <- list(
    dissimilarity_method = dissimilarity_method,
    int_val_metric = int_val_metric,
    tool_type = tool_type,
    tool_name = tool_name,
    language = language
  )

  benchmark_summary <- resource_run_tool_benchmark(prepared = prepared, tool_row = tool_row, params = params)

  # Legacy naming: external rows use "external" as dissimilarity_method and
  # the tool name as consistency_metric, matching ISC_benchmark's own scheme.
  consistency_metric_legacy <- if (identical(tool_type, "internal")) int_val_metric else tool_name
  dissimilarity_method_legacy <- if (identical(tool_type, "internal")) dissimilarity_method else "external"

  output <- data.frame(
    duration_ms = benchmark_summary$duration_ms,
    peak_memory_MB = benchmark_summary$peak_memory_MB,
    duration = benchmark_summary$duration_ms / 1000,
    memory_usage_MB = benchmark_summary$peak_memory_MB,
    cpu_usage = benchmark_summary$cpu_usage,
    tool_type = tool_type,
    tool_name = tool_name,
    language = language,
    method = consistency_metric_legacy,
    consistency_metric = consistency_metric_legacy,
    consistency.metric = consistency_metric_legacy,
    dissimilarity_method = dissimilarity_method_legacy,
    dissimilarity.method = dissimilarity_method_legacy,
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
    tool_name = tool_name
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

  expected_n <- nrow(build_resource_tool_grid(params))
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