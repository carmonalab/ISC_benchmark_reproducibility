## Installtion
#https://github.com/igrabski/sc-SHC
#devtools::install_github("igrabski/sc-SHC")

# Returns a unique value per annotation, cannot return one for each cell type
run_scSHC <- function(scTypeEval,
                      ident = NULL,
                      var.genes = NULL,
                      parallel = FALSE,
                      cores = 1){
   
   if(class(scTypeEval) != "scTypeEval"){
      stop("input should be an scTypeEval object created with `create_scTypeEval()`")
   }
   
   if(is.null(ident)){
      ident <- scTypeEval@data$`single-cell`@ident[[1]]
   }
   
   if(is.null(var.genes)){
      var.genes <- scTypeEval@gene_lists$HVG
   }
   
   data <- scTypeEval:::get_filtered_raw_matrix(scTypeEval)$counts
   
   
   eval_cluster <- scSHC::testClusters(data = data,
                                       cluster_ids = as.character(ident),
                                       batch = NULL,
                                       var.genes = var.genes,
                                       alpha = 0.05,
                                       parallel = parallel,
                                       cores = cores,
                                       num_PCs = 30
   )
   
   ret <- mclust::adjustedRandIndex(eval_cluster[[1]], as.character(ident))
   
   return(ret)
}
