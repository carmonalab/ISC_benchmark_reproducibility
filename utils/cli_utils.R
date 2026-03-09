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

message_time <- function(m) {
  time <- format(Sys.time(), "%H:%M:%S %d-%m-%Y")
  cat("\n\n####################\n", m, " -- ", time,
      "\n####################\n\n", sep = "")
}
