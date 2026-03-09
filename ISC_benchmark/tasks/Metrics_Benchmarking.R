#!/usr/bin/env Rscript

suppressPackageStartupMessages({
   library(scTypeEval)
   library(dplyr)
   library(tidyr)
   library(ggplot2)
   library(yaml)
   library(SeuratObject)
   library(Matrix)
})

source("utils/cli_utils.R")

`%||%` <- function(a, b) if (!is.null(a)) a else b

# Load scTypeEval benchmarking wrappers (installed via inst/benchmarking)
source(system.file("benchmarking", "assays_utils.R", package = "scTypeEval"))
source(system.file("benchmarking", "Metrics_benchmarking.R", package = "scTypeEval"))



# load parameters
args <- commandArgs(trailingOnly = TRUE)
args <- parse_args(args)
if (length(args) == 0) stop("Missing command-line arguments")
dataset_path <- args[["dataset"]]
config <- yaml::read_yaml(args[["params"]])
config <- convert_lists(config)
ncores <- args[["ncores"]]
config$common$ncores <- as.numeric(ncores)
# parse idents safely
idents_arg <- args[["idents"]]
idents <- strsplit(idents_arg, ",")[[1]]
dataset.name <- tools::file_path_sans_ext(basename(dataset_path))

# Ensure output directory exists
# create timestamp dir
message_time("Creating output directory")
time <- Sys.time()
timestamp <- format(time, "%H%M%S_%d%m%Y")

output_base <- args[["out"]] %||% file.path("results", "isc", "tasks", "output")
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
if(is.null(config$common$black.list)){
   data(default_black_list, package = "scTypeEval")
   bl <- list(black.list$TCR,
              black.list$Immunoglobulins,
              black.list$Ygenes) |> unlist()
   config$common$black.list <- bl
} else {
   bl <- read.table(config$common$black.list)$V1
   config$common$black.list <- bl
}


#########################
# Run loop for all params


for(ident in idents){
   message_time(paste("Computing ", ident))
   # remove NA from ident
   sub_meta = object[["metadata"]] %>% 
      filter(!is.na(.data[[ident]]))
   sub_count_matrix = object[["counts"]][,rownames(sub_meta)]
   
   # create dir
   output_dir <- file.path(task_dir,
                           ident)
   dir.create(output_dir,
              showWarnings = F,
              recursive = T)
   # preflight params for all tasks
   preparams <- list(count_matrix = sub_count_matrix,
                     metadata = sub_meta,
                     ident = ident,
                     dir = output_dir)
   
   # loop for the different vars
   if(config$run$SplitCelltype){
      message_time("Running SplitCelltype")
      params <- c(preparams,
                  config$common,
                  config$wr.SplitCelltype)
      wr <- do.call(wr.splitCellType, params)
      # compute task metrics
      r <- wr.assayPlot(wr,
                                     type = "Constant",
                                     return.df = T)
      r <- r %>% 
         mutate(ident = ident,
                dataset = dataset.name,
                task = "SplitCelltype")
      saveRDS(r, file.path(output_dir, "SplitCelltype_summary.rds"))
      
   }
   
   if(config$run$missclassify){
      message_time("Running Missclassify")
      params <- c(preparams,
                  config$common,
                  config$wr.missclasify)
      wr <- do.call(wr.missclasify, params)
      # compute task metrics
      r <- wr.assayPlot(wr,
                                     type = "monotonic",
                                     return.df = T)
      r <- r %>% 
         mutate(ident = ident,
                dataset = dataset.name,
                task = "Missclassification")
      saveRDS(r, file.path(output_dir, "Missclassification_summary.rds"))
      
   }
   
   if(config$run$Nsamples){
      message_time("Running NSamples")
      params <- c(preparams,
                  config$common,
                  config$wr.NSamples)
      wr <- do.call(wr.NSamples, params)
      # compute task metrics
      r <- wr.assayPlot(wr,
                                     type = "Constant",
                                     return.df = T)
      r <- r %>% 
         mutate(ident = ident,
                dataset = dataset.name,
                task = "NSamples")
      saveRDS(r, file.path(output_dir, "Nsamples_summary.rds"))
   }
   
   if(config$run$NCell){
      message_time("Running NCell")
      params <- c(preparams,
                  config$common,
                  config$wr.NCell)
      wr <- do.call(wr.NCell, params)
      # compute task metrics
      r <- wr.assayPlot(wr,
                                     type = "Constant",
                                     return.df = T)
      r <- r %>% 
         mutate(ident = ident,
                dataset = dataset.name,
                task = "NCell")
      saveRDS(r, file.path(output_dir, "NCell_summary.rds"))
   }
   
   
   
   if(config$run$mergeCT){
      message_time("Running mergeCT")
      
      params <- c(preparams,
                  config$common,
                  config$wr.mergeCT)
      wr <- do.call(wr.mergeCT, params)
      # compute task metrics
      r <- wr.assayPlot(wr,
                                     type = "Constant",
                                     group = "ident",
                                     return.df = T)
      r <- r %>% 
         mutate(ident = ident,
                dataset = dataset.name,
                task = "mergeCT")
      saveRDS(r, file.path(output_dir, "mergeCT_summary.rds"))
   }
   
   if(config$run$Nct){
      message_time("Running Nct")
      params <- c(preparams,
                  config$common,
                  config$wr.Nct)
      wr <- do.call(wr.Nct, params)
      # compute task metrics
      r <- wr.assayPlot(wr,
                                     type = "Constant",
                                     return.df = T)
      r <- r %>% 
         mutate(ident = ident,
                dataset = dataset.name,
                task = "Nct")
      saveRDS(r, file.path(output_dir, "Nct_summary.rds"))
   }
   
}


message_time("Capturing session info")
if (requireNamespace("renv", quietly = TRUE)) {
   renv::snapshot(lockfile = file.path(task_dir, "renv.lock"),
                  prompt = F,
                  type = "all")
}
writeLines(capture.output(sessionInfo()),
           file.path(task_dir, "sessionInfo.txt"))

message_time("Pipeline completed!")
