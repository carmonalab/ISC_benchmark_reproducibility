#!/usr/bin/env Rscript

suppressPackageStartupMessages({
   library(scTypeEval)
   library(dplyr)
   library(tidyr)
   library(ggplot2)
   library(yaml)
   library(SeuratObject)
   library(Matrix)
   library(BiocParallel)
   library(bench)
})

source("utils/cli_utils.R")
source("utils/benchmark_utils.R")

`%||%` <- function(a, b) if (!is.null(a)) a else b



# load parameters
args <- commandArgs(trailingOnly = TRUE)
args <- parse_args(args)
if (length(args) == 0) stop("Missing command-line arguments")
dataset_path <- args[["dataset"]]
config <- yaml::read_yaml(args[["params"]])
config <- convert_lists(config)
ncores <- args[["ncores"]]
# parse idents safely
idents_arg <- args[["idents"]]
idents <- strsplit(idents_arg, ",")[[1]]
dataset.name <- tools::file_path_sans_ext(basename(dataset_path))

# Ensure output directory exists
# create timestamp dir
message_time("Creating output directory")
time <- Sys.time()
timestamp <- format(time, "%H%M%S_%d%m%Y")

output_base <- args[["out"]] %||% file.path("results", "isc", "resources", "output")
if (!dir.exists(output_base)) dir.create(output_base, recursive = TRUE, showWarnings = FALSE)

task_dir <- file.path(output_base, paste(dataset.name, timestamp, sep = "_"))
if (!dir.exists(task_dir)) dir.create(task_dir, recursive = TRUE, showWarnings = FALSE)
message_time(paste("Created output directory", task_dir))

message_time("Saving parameters")
file.copy(from = args[["params"]], 
          file.path(task_dir, "params.yaml")
)


message_time("Loading datasets")
# load object and create count matrix and metadata
object <- scTypeEval::load_singleCell_object(dataset_path)


# get black list if default
if(is.null(config$black.list)){
   data(default_black_list, package = "scTypeEval")
   bl <- list(black.list$TCR,
              black.list$Immunoglobulins,
              black.list$Ygenes) |> unlist()
} else {
   bl <- read.table(config$black.list)$V1
}


combis <- expand.grid(config$dissimilarity.method,
                      config$IntVal.metric,
                      stringsAsFactors = F)

param <- BiocParallel::MulticoreParam(workers = ncores,
                                      progressbar = FALSE)

#########################
# Run loop for all params


for(ident in idents){
   message_time(paste("Computing ", ident))
   # remove NA from ident
   sub_meta = object[["metadata"]] %>% 
      filter(!is.na(.data[[ident]]))
   sub_count_matrix = object[["counts"]][,rownames(sub_meta)]
   
   #create sceval object
   sc <- create.scTypeEval(matrix = sub_count_matrix,
                           metadata = sub_meta,
                           black.list = bl)
   
   message_time("Processing data")
   sc <- Run.ProcessingData(sc,
                            ident = ident,
                            sample = config$sample,
                            normalization.method = config$normalization.method,
                            min.samples = config$min.samples,
                            min.cells = config$min.cells,
                            verbose = config$verbose)
   
   # get the gene list, the same for every run
   if(is.null(config$gene.list)){
      message_time("Computing HVG")
      sc <- Run.HVG(sc,
                    ncores = ncores)
   } else {
      sc <- add.GeneList(sc, gene.list = gene.list)
   }
   
   message_time("Computing PCA")
   sc <- Run.PCA(sc)
   
   nfeatures <- length(sc@gene.lists[[1]])
   gene.list <- names(sc@gene.lists[1])
   ncells <- ncol(sub_count_matrix)
   nsamples <- length(unique(sub_meta[[config$sample]]))
   
   mat0 <- sub_count_matrix[sc@gene.lists[[gene.list]], ]
   total_entries <- prod(dim(mat0))
   n_nonzero <- length(mat0@x)  # Only non-zero entries are stored
   sparsity <- (total_entries - n_nonzero) / total_entries
   rm(mat0)
   
   
   ili <- bplapply(1:nrow(combis),
                   BPPARAM = param,
                   function(i){
                      dm <- combis[i, 1]
                      im <- combis[i, 2]
                      
                      message(".   Running ", dm, " and ", im, "     .\n")
                      
                      rs <- benchmark_wrapper(
                         expr = {
                            
                            sc.tmp <- Run.Dissimilarity(scTypeEval = sc,
                                                        method = dm,
                                                        reduction = config$reduction,
                                                        gene.list = NULL,
                                                        black.list = NULL,
                                                        BestHit.classifier = "SingleR",
                                                        ncores = 1,
                                                        bparam = NULL,
                                                        progressbar = FALSE,
                                                        verbose = config$verbose
                            )
                            
                            sc.tmp <- get.Consistency(scTypeEval = sc.tmp,
                                                      dissimilarity.slot = dm,
                                                      Consistency.metric = im,
                                                      KNNGraph_k = config$KNNGraph_k,
                                                      hclust.method = config$hclust.method,
                                                      normalize = F,
                                                      verbose = config$verbose)
                         }
                      )
                      
                      r <- data.frame(duration = rs$duration,
                                      memory_usage_MB = rs$memory / 1024^2,
                                      cpu_usage = rs$cpu_usage,
                                      method = im,
                                      dissimilarity.method = dm)
                      return(r)
                   })
   
   dfi <- dplyr::bind_rows(ili)
   
   dfi <- dfi %>% 
      mutate(dataset = dataset.name,
             ident = ident,
             nfeatures = nfeatures,
             ncells = ncells,
             nsamples = nsamples,
             sparsity = sparsity,
             gene.list = gene.list)
   
   message_time("Saving output file")
   saveRDS(dfi,
           file.path(task_dir,
                     paste0(dataset.name, "_",
                            ident,
                            ".rds")))
}



message_time("Capturing session info")
if (requireNamespace("renv", quietly = TRUE)) {
   renv::snapshot(lockfile = file.path(task_dir,
                                       "renv.lock"),
                  prompt = F,
                  type = "all")
}
writeLines(capture.output(sessionInfo()),
           file.path(task_dir,
                     "sessionInfo.txt"))

message_time("Pipeline completed!")




