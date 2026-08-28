library(anndataR)
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
    file.path(repo_root, ".venv_311", "bin", "python"),
    file.path(repo_root, ".venv_anticor_feat", "bin", "python"),
    file.path(repo_root, ".venv_cellhint", "bin", "python")
  )

  for (candidate in candidates) {
    if (file.exists(candidate)) {
      return(candidate)
    }
  }

  reticulate_python <- Sys.getenv("RETICULATE_PYTHON", unset = "")
  if (nzchar(reticulate_python)) {
    return(reticulate_python)
  }

  "python"
}

.default_sccaf_python_bin <- function(repo_root) {
  candidate <- file.path(repo_root, ".venv_sccaf", "bin", "python")
  if (file.exists(candidate)) {
    return(candidate)
  }

  .default_python_bin(repo_root)
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


.append_flag <- function(args, flag, value) {
  if (is.null(value)) {
    return(args)
  }
  c(args, flag, as.character(value))
}


.append_bool_flag <- function(args, flag, value) {
  if (isTRUE(value)) {
    return(c(args, flag))
  }
  args
}


pipeline_spec_sccaf <- function(
  name = "sccaf",
  cluster_key = NULL,
  n = 100,
  extra_args = character()
) {
  args <- character()
  args <- .append_flag(args, "--cluster-key", cluster_key)
  args <- .append_flag(args, "--n", n)
  args <- c(args, extra_args)

  list(name = name, script = "sccaf.py", args = args)
}


pipeline_spec_cellhint <- function(
  name = "cellhint",
  cluster_key = NULL,
  dataset_key = NULL,
  random_state = 2,
  extra_args = character()
) {
  args <- character()
  args <- .append_flag(args, "--cluster-key", cluster_key)
  args <- .append_flag(args, "--dataset-key", dataset_key)
  args <- .append_flag(args, "--random-state", random_state)
  args <- c(args, extra_args)

  list(name = name, script = "cellhint.py", args = args)
}


pipeline_spec_anticor_features <- function(
  name = "anticor_features",
  cluster_key = NULL,
  min_cells = 10,
  max_cell_types = NULL,
  species = "hsapiens",
  n_rand_feat = 2000,
  fpr = 0.001,
  fdr = 1 / 15,
  num_pos_cor = 10,
  bin_size = 5000,
  scratch_dir = NULL,
  offline_mode = FALSE,
  use_live_pathway_lookup = FALSE,
  score_k = 1.0,
  extra_args = character()
) {
  args <- character()
  args <- .append_flag(args, "--cluster-key", cluster_key)
  args <- .append_flag(args, "--min-cells", min_cells)
  args <- .append_flag(args, "--max-cell-types", max_cell_types)
  args <- .append_flag(args, "--species", species)
  args <- .append_flag(args, "--n-rand-feat", n_rand_feat)
  args <- .append_flag(args, "--fpr", fpr)
  args <- .append_flag(args, "--fdr", fdr)
  args <- .append_flag(args, "--num-pos-cor", num_pos_cor)
  args <- .append_flag(args, "--bin-size", bin_size)
  args <- .append_flag(args, "--scratch-dir", scratch_dir)
  args <- .append_bool_flag(args, "--offline-mode", offline_mode)
  args <- .append_bool_flag(args, "--use-live-pathway-lookup", use_live_pathway_lookup)
  args <- .append_flag(args, "--score-k", score_k)
  args <- c(args, extra_args)

  list(name = name, script = "anticor_features.py", args = args)
}


pipeline_spec_popv <- function(
  name = "popv",
  reference = NULL,
  ref_labels_key = NULL,
  ref_batch_key = NULL,
  query_batch_key = NULL,
  query_celltype_key = NULL,
  prediction_mode = "fast",
  methods = NULL,
  pretrained_model_repo = "popV/tabula_sapiens_All_Cells",
  hub_cache_dir = NULL,
  n_samples_per_label = 300,
  hvg = 4000,
  unknown_celltype_label = "unknown",
  save_path_trained_models = "tmp/popv_models",
  extra_args = character()
) {
  args <- character()
  args <- .append_flag(args, "--reference", reference)
  args <- .append_flag(args, "--ref-labels-key", ref_labels_key)
  args <- .append_flag(args, "--ref-batch-key", ref_batch_key)
  args <- .append_flag(args, "--query-batch-key", query_batch_key)
  args <- .append_flag(args, "--query-celltype-key", query_celltype_key)
  args <- .append_flag(args, "--prediction-mode", prediction_mode)
  if (!is.null(methods)) {
    method_value <- if (length(methods) > 1) paste(methods, collapse = ",") else methods
    args <- .append_flag(args, "--methods", method_value)
  }
  args <- .append_flag(args, "--pretrained-model-repo", pretrained_model_repo)
  args <- .append_flag(args, "--hub-cache-dir", hub_cache_dir)
  args <- .append_flag(args, "--n-samples-per-label", n_samples_per_label)
  args <- .append_flag(args, "--hvg", hvg)
  args <- .append_flag(args, "--unknown-celltype-label", unknown_celltype_label)
  args <- .append_flag(args, "--save-path-trained-models", save_path_trained_models)
  args <- c(args, extra_args)

  list(name = name, script = "popv.py", args = args)
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

  adata <- anndataR::AnnData(
    X = Matrix::t(filt_data$counts),
    obs = as.data.frame(filt_data$metadata)
  )
  anndataR::write_h5ad(adata, h5ad_path)
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
  required_output_cols = c("celltype", "score"),
  continue_on_error = FALSE,
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
  pipeline_errors <- list()
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

    pipeline_python_bin <- python_bin
    if (identical(spec$name, "sccaf")) {
      pipeline_python_bin <- .default_sccaf_python_bin(repo_root)
    }

    if (!file.exists(pipeline_python_bin)) {
      stop(
        "Python interpreter not found for ", spec$name, ": ", pipeline_python_bin
      )
    }
    if (file.access(pipeline_python_bin, mode = 1) != 0) {
      stop(
        "Python interpreter is not executable for ", spec$name, ": ", pipeline_python_bin
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

    cmd_display <- paste(
      c(shQuote(pipeline_python_bin), shQuote(cmd_args)),
      collapse = " "
    )

    stdout_log <- tempfile(pattern = "py_pipeline_stdout_", fileext = ".log")
    stderr_log <- tempfile(pattern = "py_pipeline_stderr_", fileext = ".log")

    run_status <- tryCatch(
      suppressWarnings(
        system2(
          command = pipeline_python_bin,
          args = cmd_args,
          stdout = stdout_log,
          stderr = stderr_log
        )
      ),
      error = function(e) e
    )

    read_log_file <- function(path) {
      if (!file.exists(path)) {
        return(character())
      }
      lines <- readLines(path, warn = FALSE)
      if (length(lines) == 0) {
        return(character())
      }
      lines
    }

    stdout_lines <- read_log_file(stdout_log)
    stderr_lines <- read_log_file(stderr_log)
    unlink(c(stdout_log, stderr_log), force = TRUE)

    if (inherits(run_status, "error")) {
      failure_message <- paste0(
        "Python pipeline failed for ", spec$name, " (", basename(script_path), ").\n",
        "Launcher error: ", conditionMessage(run_status), "\n",
        "Command: ", cmd_display, "\n",
        "stdout:\n", paste(stdout_lines, collapse = "\n"), "\n",
        "stderr:\n", paste(stderr_lines, collapse = "\n")
      )
      if (isTRUE(continue_on_error)) {
        message("[external] ", failure_message)
        pipeline_errors[[spec$name]] <- failure_message
        next
      }
      stop(failure_message)
    }

    exit_code <- as.integer(run_status)

    if (!identical(exit_code, 0L)) {
      signal_note <- ""
      if (!is.na(exit_code) && exit_code >= 128L) {
        signal_note <- paste0(" (possible signal ", exit_code - 128L, ")")
      }
      failure_message <- paste0(
        "Python pipeline failed for ", spec$name, " (", basename(script_path), ").\n",
        "Exit code: ", exit_code, signal_note, "\n",
        "Command: ", cmd_display, "\n",
        "stdout:\n", paste(stdout_lines, collapse = "\n"), "\n",
        "stderr:\n", paste(stderr_lines, collapse = "\n")
      )
      if (isTRUE(continue_on_error)) {
        message("[external] ", failure_message)
        pipeline_errors[[spec$name]] <- failure_message
        next
      }
      stop(failure_message)
    }

    if (!file.exists(output_csv)) {
      failure_message <- paste0(
        "Pipeline completed but output CSV was not created for ", spec$name, ": ", output_csv
      )
      if (isTRUE(continue_on_error)) {
        message("[external] ", failure_message)
        pipeline_errors[[spec$name]] <- failure_message
        next
      }
      stop(failure_message)
    }

    parsed_result <- tryCatch(
      read_csv(
        output_csv,
        stringsAsFactors = FALSE,
        check.names = FALSE
      ),
      error = function(e) e
    )
    if (inherits(parsed_result, "error")) {
      failure_message <- paste0(
        "Pipeline output parsing failed for ", spec$name, ": ", conditionMessage(parsed_result)
      )
      if (isTRUE(continue_on_error)) {
        message("[external] ", failure_message)
        pipeline_errors[[spec$name]] <- failure_message
        next
      }
      stop(failure_message)
    }
    results[[spec$name]] <- parsed_result

    missing_cols <- setdiff(required_output_cols, colnames(results[[spec$name]]))
    if (length(missing_cols) > 0) {
      failure_message <- paste0(
        "Pipeline output is missing required column(s) for ", spec$name, ": ",
        paste(missing_cols, collapse = ", "),
        ". Expected at least: ", paste(required_output_cols, collapse = ", "),
        ". Output file: ", output_csv
      )
      if (isTRUE(continue_on_error)) {
        message("[external] ", failure_message)
        pipeline_errors[[spec$name]] <- failure_message
        next
      }
      stop(failure_message)
    }

    result_paths[[spec$name]] <- output_csv

    if (isTRUE(assign_to_env)) {
      assign(spec$name, results[[spec$name]], envir = envir)
    }
  }

  attr(results, "input_h5ad") <- h5ad_path
  attr(results, "output_csvs") <- result_paths
  attr(results, "pipeline_errors") <- pipeline_errors
  results
}
