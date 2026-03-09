# Helpers to make all scripts runnable from the repo root.

project_root <- function() {
  # When sourcing, getwd() should already be the project root (recommended).
  # This helper exists mostly to centralize conventions.
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

path_from_root <- function(...) {
  file.path(project_root(), ...)
}

ensure_dir <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
  invisible(path)
}
