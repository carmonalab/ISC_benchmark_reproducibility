convert_lists <- function(x) {
  if (is.list(x)) {
    x <- lapply(x, convert_lists)
    if (is.null(names(x)) && all(sapply(x, is.atomic))) {
      return(unlist(x, use.names = FALSE))
    }
  }
  x
}

parse_args <- function(args) {
  arg_list <- list()
  if (length(args) %% 2 != 0) {
    stop("Invalid CLI args: expected --key value pairs")
  }
  for (i in seq(1, length(args), 2)) {
    key <- gsub("^--", "", args[[i]])
    val <- args[[i + 1]]
    arg_list[[key]] <- val
  }
  arg_list
}

# ============================================================================
# MESSAGE HELPERS
# ============================================================================

message_time <- function(...) {
  msg <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", ...)
  message(msg)
}

message_step <- function(step_name, ...) {
  msg <- paste0("[", get_current_pipeline(), "] ", step_name, " — ", ...)
  message_time(msg)
}

# ============================================================================
# METADATA HELPERS
# ============================================================================

#' Normalize metadata column names by removing special characters
#'
#' Standardizes string values by:
#' 1. Replacing non-alphanumeric characters (except . and -) with dots
#' 2. Collapsing consecutive dots into single dots
#' 3. Trimming whitespace from both ends
#' 4. Removing trailing dots
#'
#' Useful for cleaning batch, condition, or tissue annotations to ensure
#' valid filenames and consistent naming conventions.
#'
#' @param x Character vector to normalize
#'
#' @return Character vector with normalized values
#'
#' @examples
#' \dontrun{
#' normalize_metadata_name(c("tissue sample", "pre-post", "2D/3D"))
#' # Returns: c("tissue.sample", "pre-post", "2D.3D")
#' }
normalize_metadata_name <- function(x) {
  if (is.factor(x)) x <- as.character(x)
  if (!is.character(x)) return(x)
  
  # Replace non-alphanumeric (except . and -) with dots
  x <- stringr::str_replace_all(x, "[^a-zA-Z0-9._-]", ".")
  # Collapse consecutive dots
  x <- stringr::str_replace_all(x, "\\.+", ".")
  # Trim whitespace
  x <- stringr::str_trim(x, side = "both")
  # Remove trailing dots
  x <- sub("\\.$", "", x)
  
  x
}
