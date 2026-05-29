#' Perturb counts for selected cells in a raw counts matrix
#'
#' @param counts A raw counts matrix (genes x cells)
#' @param colnames A character vector of cell (column) names to perturb
#' @param noise Relative noise level (e.g., 0.05 = low, 0.6 = high)
#' @return A modified counts matrix with perturbed cells
#' @examples
#' perturbed <- perturb_cells(counts, c("cell1", "cell2"), noise = 0.3)
perturb_cells <- function(counts, colnames, noise = 0.5) {
  stopifnot(is.matrix(counts) || is.data.frame(counts))
  stopifnot(all(colnames %in% colnames(counts)))
  stopifnot(is.numeric(noise) && length(noise) == 1 && noise >= 0)
  counts <- as.matrix(counts)
  set.seed(42) # for reproducibility
  for (cell in colnames) {
    # Add random noise proportional to counts
    original <- counts[, cell]
    # Poisson noise: add/subtract random counts
    perturb <- rpois(length(original), lambda = original * noise)
    sign <- sample(c(-1, 1), length(original), replace = TRUE)
    new_counts <- original + sign * perturb
    new_counts[new_counts < 0] <- 0
    counts[, cell] <- new_counts
  }
  return(counts)
}
