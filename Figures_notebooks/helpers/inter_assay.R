get_method_sample_combinations <- function(metadata,
                                           sample_id = "Sample",
                                           batch_id = "method",
                                           seed = 22){
   rlist <- list()
   # get unique sample-batch
   umd <- metadata %>%
      distinct(.data[[sample_id]], .data[[batch_id]])
   #umd[[sample_id]] <- scTypeEval:::purge_label(umd[[sample_id]])
   # all samples for all batches
   rlist[["all"]] <- umd[[sample_id]]
   # all samples from each batch
   spl <- split(umd[[sample_id]], umd[[batch_id]])
   
   
   set.seed(seed)
   spl_set <- lapply(spl, function(x){
      l <- floor(length(x)/2)
      a <- sample(x, size = l)
      b <- x[!x %in% a]
      r <- list(A = a,
                B =b)
   })
   
   rlist <- c(rlist, spl,
              unlist(spl_set, recursive = F))
   
   return(rlist)
      
}

repo_root <- if (dir.exists("sample_agnostic_utils")) {
   getwd()
} else {
   normalizePath("..", mustWork = TRUE)
}
source(file.path(repo_root, "sample_agnostic_utils", "run_python_pipelines.R"))


run_isc <- function(counts,
                    metadata,
                    which_samples = NULL,
                    ident = "cell_type",
                    sample = "Sample",
                    min_samples = 3,
                    min_cells = 10,
                    dm = c("Pseudobulk:Cosine",
                           "recip_classif:Match"),
                    black_list = NULL,
                    hvg = NULL
) {
   meta <- metadata %>% 
      filter(.data[[sample]] %in% which_samples)
   cc <- counts[,rownames(meta)]
   
   sc <- create_scTypeEval(matrix = cc,
                           metadata = meta,
                           black_list = bl)
   message("> Running processing\n")
   sc <- run_processing_data(sc,
                             ident = ident,
                             sample = sample,
                             # enforce more cell type by lower threshold of samples
                             min_samples = min_samples,
                             min_cells = min_cells,
                             verbose = F)
   if(is.null(hvg)){
      sc <- run_hvg(sc, 
                    aggregation = "single-cell",
                    ngenes = 2000,
                    var_method = "scran",
                    ncores = 5, 
                    verbose = F)
   } else {
      sc <- add_gene_list(sc, gene_list = list(hvg = hvg))
   }
   
   # no need to run pca on single-cell data
   single_mat <- sc@data$`single-cell`
   sc@data$`single-cell` <- NULL
   message("> Running PCA\n")
   sc <- run_pca(sc,
                 verbose = F)
   
   for(mdiss in dm){
      cat("#### Running ", mdiss, " ####\n")
      sc <- run_dissimilarity(sc,
                              method = mdiss,
                              ncores = 7,
                              verbose = F)
   }
   # return single-cell data for SCCAF and MetaNeighbor
   sc@data$`single-cell` <- single_mat
   
   message("> Running consistency\n")
   cons <- get_consistency(sc,
                           consistency_metric = c("silhouette",
                                                  "2label_silhouette",
                                                  "MetaNeighbor_Supervised")
   )
   
   # Run SCCAF
   message("> Running SCCAF\n")
   sccaf_results <- run_sample_agnostic_python_pipelines(
      scTypeEval = sc,
      pipelines = list(
         sccaf = pipeline_spec_sccaf(cluster_key = ident, n = 100)
      ),
      file_prefix = basename(tempfile("thymus_")),
      continue_on_error = FALSE,
      cleanup = TRUE
   )$sccaf
   
   
   sccaf <- data.frame(celltype = sccaf_results$celltype,
                       measure = sccaf_results$f1_score,
                       consistency_metric = "SCCAF",
                       dissimilarity_method = NA,
                       ident = sccaf_results$cluster_key)
   cons <- bind_rows(cons, sccaf)
   
   # make scTypeEval object lighter
   sc@data$`single-cell` <- NULL
   sc@counts <- Matrix::sparseMatrix(
      i = integer(0),
      j = integer(0),
      x = numeric(0),
      dims = dim(sc@counts),
      dimnames = dimnames(sc@counts)
   )
   sc@metadata <- data.frame()
   
   
   return(list(scTypeEval = sc,
               consistency = cons))
}
