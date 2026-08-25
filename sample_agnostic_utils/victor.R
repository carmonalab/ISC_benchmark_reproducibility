# Installation
# https://github.com/Charlene717/VICTOR
# devtools::install_github("Charlene717/VICTOR")

library(Seurat)
library(VICTOR)
library(harmony)


get_seurat <- function(scTypeEval, prefix = "ident") {
   filt <- scTypeEval:::get_filtered_raw_matrix(scTypeEval)
   
   ident_name <- names(scTypeEval@data$`single-cell`@ident)
   md <- filt$metadata[, ident_name, drop = FALSE]
   names(md) <- prefix
   seu <- SeuratObject::CreateSeuratObject(
      counts = filt$counts,
      meta.data = md
   )
   seu <- Seurat::NormalizeData(
      seu,
      verbose = FALSE
   )
   seu <- Seurat::FindVariableFeatures(
      seu,
      selection.method = "vst",
      nfeatures = min(2000, nrow(seu)),
      verbose = FALSE
   )
   seu <- Seurat::ScaleData(
      seu,
      features = Seurat::VariableFeatures(seu),
      verbose = FALSE
   )
   seu <- Seurat::RunPCA(
      seu,
      features = Seurat::VariableFeatures(seu),
      npcs = 50,
      verbose = FALSE
   )
   return(seu)
}

run_victor <- function(query,
                       ref){
   seuratObject_Query <- get_seurat(query, prefix = "Annotation")
   seuratObject_Ref <- get_seurat(ref, prefix = "Actual_Cell_Type")
   
   
   lt <- VICTOR::VICTOR(seuratObject_Query,
                        seuratObject_Ref,
                        ActualCellTypeColumn = "Actual_Cell_Type",
                        AnnotCellTypeColumn = "Annotation")
}
