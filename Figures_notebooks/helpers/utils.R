reverseNormalize <- function(cell) {
   # reverse logn and substract 1
   cellE <- expm1(cell)
   #value for 1, most common no 0 value
   sc <- min(cellE[cellE>0])
   s <- cellE / sc
   # round because we lose some distant decimals in the process
   s <- as.integer(round(s))
   names(s) <- names(cellE)
   return(s)
}
reverseNormalize_matrix <- function(x) {
   x <- methods::as(x, "dgCMatrix")
   
   for (j in seq_len(ncol(x))) {
      # Identify the positions in x@x corresponding to nonzero values
      # in column j. The x@p slot stores column boundaries.
      idx <- seq.int(x@p[j] + 1L, x@p[j + 1L])
      # Skip cells with no stored nonzero values.
      if (length(idx) == 0L) next
      # Reverse the log1p transformation: log(1 + value) -> value.
      v <- expm1(x@x[idx])
      # Identify positive values to calculate the cell-specific scaling factor.
      pos <- v > 0
      
      if (any(pos)) {
         # Use the smallest positive value as the scaling factor.
         # Divide by it and round to obtain integer-like values.
         v <- round(v / min(v[pos]))
         # Replace the stored expression values in the sparse matrix.
         x@x[idx] <- v
      } else {
         # If no positive values remain, set the stored values to zero.
         x@x[idx] <- 0
      }
   }
   
   Matrix::drop0(x)
}


check_counts <- function(x) {
   vals <- if (inherits(x, "sparseMatrix")) x@x else as.vector(x)
   vals <- vals[is.finite(vals)]
   
   cat("Integer-valued:", all(abs(vals - round(vals)) < 1e-8), "\n")
   cat("Range:", range(vals), "\n")
   cat("Quantiles:\n")
   print(quantile(vals, c(0, 0.5, 0.9, 0.99, 1)))
   
   invisible(NULL)
}


downsample_df <- function(df, col1, col2, n = 200, seed = 22) {
   set.seed(seed)
   
   out <- df %>%
      tibble::rownames_to_column("cellid") %>% 
      dplyr::group_by(.data[[col1]], .data[[col2]]) %>%
      dplyr::group_modify(~ dplyr::slice_sample(
         .x,
         n = min(nrow(.x), n),
         replace = FALSE
      )) %>%
      dplyr::ungroup() %>% 
      tibble::column_to_rownames("cellid")
   
   
   return(out)
}
