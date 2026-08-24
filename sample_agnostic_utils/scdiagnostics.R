# Installation
# https://github.com/ccb-hms/scDiagnostics
# !! Requires R >4.4

get_sce_subset <- function(scTypeEval,
                           filter_genes = FALSE){
   
   
   if(filter_genes){
      hvg <- scTypeEval@gene_lists$HVG
   } else {
      hvg <- rownames(scTypeEval@data$`single-cell`@matrix)
   }
   
   logcounts <- scTypeEval@data$`single-cell`@matrix[hvg,]
   cell_metadata <- data.frame(rownames = colnames(logcounts),
                               ident = as.character(scTypeEval@data$`single-cell`@ident[[1]]))
   
   sce <- SingleCellExperiment::SingleCellExperiment(
      assays = list(logcounts = logcounts),
      colData = cell_metadata
      )
   return(sce)
}


run_scdiagnostics <- function(query,
                              ref){
   

   reference_data <- get_sce_subset(scTypeEval = ref)
   reference_data <- scDiagnostics::processPCA(reference_data)
   query_data <- get_sce_subset(scTypeEval = query)
   query_data <- scDiagnostics::processPCA(query_data)
   
   cramer <- scDiagnostics::calculateCramerPValue(
      reference_data = reference_data,
      query_data = query_data,
      ref_cell_type_col = "ident",
      query_cell_type_col = "ident",
      cell_types = NULL, # ALL
      assay_name = "logcounts",
      pc_subset = 1:5,
      max_cells_ref = 5000,
      max_cells_query = 5000
   )
   
   hotelling <- scDiagnostics::calculateHotellingPValue(
      reference_data = reference_data,
      query_data = query_data,
      ref_cell_type_col = "ident",
      query_cell_type_col = "ident",
      cell_types = NULL, # ALL
      pc_subset = 1:5,
      n_permutation = 500,
      assay_name = "logcounts",
      max_cells_query = 5000,
      max_cells_ref = 5000
   )
   
   ret <- data.frame(cramer = 1- cramer,
                     hotelling = 1 - hotelling)
   ret$celltype <- rownames(ret)
   
   return(ret)
   
}
