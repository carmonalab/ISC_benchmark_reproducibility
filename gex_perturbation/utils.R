#' Perturb counts for selected cells in a raw counts matrix
#'
#' @param counts A raw counts matrix (genes x cells)
#' @param colnames A character vector of cell (column) names to perturb
#' @param noise Relative noise level (e.g., 0.05 = low, 0.6 = high)
#' @return A modified counts matrix with perturbed cells
#' @examples
#' perturbed <- perturb_cells(counts, c("cell1", "cell2"), noise = 0.3)
perturb_cells <- function(counts, colnames, noise = 0.5, seed = 22) {
  stopifnot(all(colnames %in% colnames(counts)))
  stopifnot(is.numeric(noise) && length(noise) == 1 && noise >= 0)
  counts <- as.matrix(counts)
  set.seed(seed) # for reproducibility
  # Vectorized perturbation using apply
  counts[, colnames] <- apply(counts[, colnames, drop = FALSE], 2, function(original) {
    perturb <- rpois(length(original), lambda = original * noise)
    sign <- sample(c(-1, 1), length(original), replace = TRUE)
    new_counts <- original + sign * perturb
    new_counts[new_counts < 0] <- 0
    return(new_counts)
  })
  return(counts)
}


get_consistency_wide <- function(scTypeEval,
                                 consistency_metric = c("silhouette", "2label_silhouette"),
                                 global = "2label_silhouette | Pseudobulk:Cosine",
                                 local = "silhouette | recip_classif:Match"){
   
   consis <- get_consistency(scTypeEval,
                             consistency_metric = consistency_metric)
   consis <- consis %>% 
      mutate(method_type = paste(consistency_metric,
                                 dissimilarity_method,
                                 sep = " | ")) %>% 
      filter(method_type %in% c(global, local)) %>% 
      mutate(method_type = factor(method_type,
                                  levels = c(global, local),
                                  labels = c("global", "local")))
   wide <- consis %>% 
      select(method_type, measure, celltype) %>% 
      pivot_wider(names_from = "method_type",
                  values_from = "measure") %>% 
      mutate(product = global * local)
   return(wide)
}
