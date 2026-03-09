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

# Plot helpers used by this script
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
split_col <- args[["split_col"]] %||% "split"
# parse idents safely
idents_arg <- args[["idents"]]
idents <- trimws(strsplit(idents_arg, ",")[[1]])
idents <- idents[idents != ""]
dataset.name <- tools::file_path_sans_ext(basename(dataset_path))

# Ensure output directory exists
# create timestamp dir
message_time("Creating output directory")
time <- Sys.time()
timestamp <- format(time, "%H%M%S_%d%m%Y")

task_dir <- paste(dataset.name, timestamp, sep = "_")
output_base <- args[["out"]] %||% file.path("results", "isc", "batch_effect_singleRun", "output")
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

if (!split_col %in% colnames(object$metadata)) {
   stop("Grouping column not found in object$metadata: ", split_col)
}


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
# generate combinations of datasets
batches <- unique(object$metadata[[split_col]])
batches_comb <- scTypeEval:::sample_variable_length_combinations(batches,
                                                                 min_k = 2,
                                                                 max_k = 2,
                                                                 num_samples = 60,
                                                                 seed = 22)
# add singles batches
#for(b in batches){batches_comb[[b]] <- b}

# loop with all datasets combinations

for(b in rev(names(batches_comb))){
   
   batch_dir <- file.path(task_dir, b)
   if (!dir.exists(batch_dir)) {
      dir.create(batch_dir, recursive = TRUE)
   }
   message_time(paste("Created output directory", batch_dir))
   
   
   message_time(paste("Subsetting object for", b, "group"))
   metadata <- object$metadata |>
      filter(.data[[split_col]] %in% batches_comb[[b]])
   matrix <- object$counts[, rownames(metadata)]
   
   
   
   # loop for each ident
   
   for(ident in idents){
      message_time(paste("Computing ", ident))
      # remove NA from ident
      sub_meta = metadata %>% 
         filter(!is.na(.data[[ident]]))
      sub_count_matrix = matrix[,rownames(sub_meta)]
      
      if(length(unique(sub_meta[[ident]]))<2){ next }
      
      # create dir
      output_dir <- file.path(batch_dir,
                              ident)
      dir.create(output_dir,
                 showWarnings = F,
                 recursive = T)
      # preflight params for all tasks
      preparams <- list(count_matrix = sub_count_matrix,
                        metadata = sub_meta,
                        ident = ident)
      
      params <- c(preparams,
                  config$common)
      
      sc <- do.call(scTypeEval:::wrapper_scTypeEval, params)
      
      res <- list(Dissimilarity = NULL,
                  Consistency = NULL)
      
      diss <- sc@dissimilarity
      diss <- lapply(diss, function(x){
         as.matrix(x@dissimilarity)
      })
      res[["Dissimilarity"]] <- diss
      
      # --- Safe call to get.Consistency (error only) ---
      params <- c(list(scTypeEval = sc),
                  config$consistency)
      cons <- tryCatch(
         do.call(get.Consistency, params) %>% 
            mutate(batch = b,
                   dataset = dataset.name),
         error = function(e) {
            message("get.Consistency() failed for ", ident, ": ", e$message)
            NULL
         }
      )
      
      res[["Consistency"]] <- cons
      
      saveRDS(res, file.path(batch_dir, paste0(b, "_", ident, ".rds")))
      wrapper_plots(sc,
                                 dir.path = file.path(batch_dir, paste0(b, "_", ident)),
                                 reduction = config$common$reduction)
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



