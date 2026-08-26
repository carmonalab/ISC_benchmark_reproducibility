library(anndata)
library(Matrix)

.find_repo_root <- function(start_dir = getwd()) {
  current <- normalizePath(start_dir, mustWork = TRUE)

  repeat {
    if (dir.exists(file.path(current, "sample_agnostic_utils"))) {
      return(current)
    }

    parent <- dirname(current)
    if (identical(parent, current)) {
      stop(
        "Could not find repository root containing sample_agnostic_utils/ from: ",
        start_dir
      )
    }

    current <- parent
  }
}


.default_python_bin <- function(repo_root) {
  candidates <- c(
    file.path(repo_root, ".venv", "bin", "python"),
    file.path(repo_root, ".venv_311", "bin", "python")
  )

  for (candidate in candidates) {
    if (file.exists(candidate)) {
      return(normalizePath(candidate, mustWork = TRUE))
    }
  }

  reticulate_python <- Sys.getenv("RETICULATE_PYTHON", unset = "")
  if (nzchar(reticulate_python)) {
    return(reticulate_python)
  }

  "python"
}


.normalize_pipeline_specs <- function(pipelines) {
  if (is.character(pipelines)) {
    pipelines <- as.list(pipelines)
  }

  if (!is.list(pipelines) || length(pipelines) == 0) {
    stop("pipelines must be a non-empty character vector or list.")
  }

  specs <- vector("list", length(pipelines))
  spec_names <- names(pipelines)

  for (i in seq_along(pipelines)) {
    spec <- pipelines[[i]]

    if (is.character(spec) && length(spec) == 1) {
      script <- spec
      args <- character()
    } else if (is.list(spec) && !is.null(spec$script)) {
      script <- spec$script
      args <- spec$args %||% character()
    } else {
      stop(
        "Each pipeline must be either a script name/path or a list with `script` and optional `args`."
      )
    }

    if (!is.character(script) || length(script) != 1) {
      stop("Each pipeline script must be a single string.")
    }
    if (!is.character(args)) {
      stop("Pipeline args must be a character vector.")
    }

    name <- spec_names[[i]]
    if (is.null(name) || identical(name, "")) {
      name <- tools::file_path_sans_ext(basename(script))
    }

    specs[[i]] <- list(
      name = name,
      script = script,
      args = args
    )
  }

  specs
}


.resolve_script_path <- function(script, repo_root) {
  if (grepl("^/", script)) {
    path <- script
  } else {
    path <- file.path(repo_root, "sample_agnostic_utils", script)
  }

  normalizePath(path, mustWork = TRUE)
}


.input_flag_for_script <- function(script_path) {
  if (identical(basename(script_path), "popv.py")) {
    return("--query")
  }
  "--input"
}


.export_sctypeeval_to_h5ad <- function(scTypeEval, h5ad_path) {
  filt_data <- scTypeEval:::get_filtered_raw_matrix(scTypeEval)

  adata <- anndata::AnnData(
    X = Matrix::t(filt_data$counts),
    obs = as.data.frame(filt_data$metadata)
  )
  anndata::write_h5ad(adata, h5ad_path)
  h5ad_path
}


`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}


run_sample_agnostic_python_pipelines <- function(
  scTypeEval,
  pipelines,
  python_bin = NULL,
  tmp_dir = NULL,
  file_prefix = NULL,
  assign_to_env = FALSE,
  envir = parent.frame(),
  read_csv = utils::read.csv,
  cleanup = FALSE
) {
  repo_root <- .find_repo_root()

  if (is.null(python_bin) || identical(python_bin, "")) {
    python_bin <- .default_python_bin(repo_root)
  }

  if (is.null(tmp_dir)) {
    tmp_dir <- file.path(repo_root, "sample_agnostic_utils", "tmp")
  }
  dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)

  specs <- .normalize_pipeline_specs(pipelines)

  if (is.null(file_prefix) || identical(file_prefix, "")) {
    file_prefix <- paste0("sctypeeval_", format(Sys.time(), "%Y%m%d_%H%M%S"))
  }

  h5ad_path <- file.path(tmp_dir, paste0(file_prefix, ".h5ad"))
  .export_sctypeeval_to_h5ad(scTypeEval, h5ad_path)

  results <- vector("list", length(specs))
  result_paths <- vector("list", length(specs))
  names(results) <- vapply(specs, `[[`, character(1), "name")
  names(result_paths) <- names(results)

  on.exit({
    if (cleanup) {
      unlink(h5ad_path)
      unlink(unlist(result_paths, use.names = FALSE))
    }
  }, add = TRUE)

  for (i in seq_along(specs)) {
    spec <- specs[[i]]
    script_path <- .resolve_script_path(spec$script, repo_root)
    output_csv <- file.path(tmp_dir, paste0(file_prefix, "_", spec$name, ".csv"))

    reserved_args <- c("--input", "--query", "--output")
    if (any(spec$args %in% reserved_args)) {
      stop(
        "Pipeline args for ", spec$name,
        " must not include --input, --query, or --output; these are managed by the wrapper."
      )
    }

    cmd_args <- c(
      script_path,
      .input_flag_for_script(script_path),
      h5ad_path,
      "--output",
      output_csv,
      spec$args
    )

    run_log <- system2(
      command = python_bin,
      args = cmd_args,
      stdout = TRUE,
      stderr = TRUE
    )
    exit_code <- attr(run_log, "status") %||% 0L

    if (!identical(exit_code, 0L)) {
      stop(
        "Python pipeline failed for ", spec$name, " (", basename(script_path), ").\n",
        paste(run_log, collapse = "\n")
      )
    }

    if (!file.exists(output_csv)) {
      stop(
        "Pipeline completed but output CSV was not created for ", spec$name, ": ", output_csv
      )
    }

    results[[spec$name]] <- read_csv(
      output_csv,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    result_paths[[spec$name]] <- output_csv

    if (isTRUE(assign_to_env)) {
      assign(spec$name, results[[spec$name]], envir = envir)
    }
  }

  attr(results, "input_h5ad") <- h5ad_path
  attr(results, "output_csvs") <- result_paths
  results
}