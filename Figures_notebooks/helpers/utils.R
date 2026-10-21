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

reverseNormalize_matrix <- function(mat, margin = 2L) {
   
   # Ensure compressed sparse column format
   mat <- methods::as(mat, "dgCMatrix")
   
   # For margin = 1, transpose so each target vector is a column
   if (margin == 1L) {
      mat <- Matrix::t(mat)
   } else if (margin != 2L) {
      stop("`margin` must be 1L (rows) or 2L (columns)")
   }
   
   # Reverse log1p transformation on stored nonzero values only
   x <- expm1(mat@x)
   
   # Identify the column index for every stored nonzero value
   cell_idx <- rep(seq_len(ncol(mat)), diff(mat@p))
   
   # Calculate minimum positive value per cell
   min_vals <- tapply(x[x > 0], cell_idx[x > 0], min)
   
   # Map scaling factors back to stored values
   scale <- min_vals[cell_idx]
   
   # Rescale and round
   valid <- !is.na(scale)
   x[valid] <- round(x[valid] / scale[valid])
   x[!valid] <- 0L
   
   # Update sparse matrix and remove explicit zeros
   mat@x <- x
   mat <- Matrix::drop0(mat)
   
   # Restore original orientation if needed
   if (margin == 1L) {
      mat <- Matrix::t(mat)
   }
   
   mat
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
