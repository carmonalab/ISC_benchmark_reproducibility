wd <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
# load packages
library(dplyr)
library(tibble)
library(BiocParallel)
library(tidyr)
library(ggplot2)
library(scTypeEval)
library(stringr)
library(purrr)
library(data.table)

source("utils/cli_utils.R")

# Function to compute method_type ranking
compute_ranking <- function(selected_dss, alldf) {
   # Weighted average score per method_type
   full_score <- alldf %>%
      mutate(dss = paste(dataset, ident, sep = "_")) %>%
      filter(dss %in% selected_dss) %>%
      group_by(method_type, assay) %>%
      summarise(score = mean(score, na.rm = TRUE), .groups = "drop") %>%  # one score per assay
      group_by(method_type) %>%
      summarise(score = weighted.mean(score, w = weight_vector[assay], na.rm = TRUE)) %>% 
      ungroup() %>% 
      mutate(score = round(score, 3)) %>%
      arrange(desc(score)) %>%
      pull(method_type)
   
   return(full_score)
}

weight_vector <- c("Missclassification" = 0.2,
                   "Similar-celltype" = 0.2,
                   "Number-samples" = 0.1,
                   "Complexity" = 0.1,
                   "Number-cells" = 0.1,
                   "Batch-Effect" = 0.1,
                   "Perturbation-Effect" = 0.1,
                   "Granularity"= 0.1)



savedir <- file.path(wd, "results/isc/saturation/output")
dir.create(savedir,
           recursive = T,
           showWarnings = F)

##############################################################################
# Load files
##############################################################################
input_dir <- file.path(wd, "results/isc/saturation/data")
type <- "ident"
input_file <- file.path(input_dir, paste0(type, ".rds"))

alldf <- readRDS(input_file)


# Max number of disjoint pairs to collect
n_pairs <- 5e3
chunk_size <- 100
max_iter <- 200
seed <- 22

dss <- unique(paste(alldf$dataset, alldf$ident, sep = "_"))
n_total <- length(dss)
bp_param <- scTypeEval:::set_parallel_params(ncores = 16, progressbar = FALSE)

all_results <- vector("list", floor(n_total / 2))
set.seed(seed)

# Helper to convert indices into bitset
as_bit <- function(indices, n_total) {
   bit <- logical(n_total)
   bit[indices] <- TRUE
   bit
}

for (k in seq_len(13)) {
   message("Processing group size: ", k)
   
   combs_k <- combn(seq_len(n_total), m = k, simplify = FALSE)
   n_combs <- length(combs_k)
   message("Found ", n_combs, " dataset combinations.")
   
   # Precompute bitsets
   if (length(combs_k) > 500) {
      bitsets <- bplapply(combs_k, as_bit, n_total = n_total, BPPARAM = bp_param)
   } else {
      bitsets <- lapply(combs_k, as_bit, n_total = n_total)
   }
   
   valid_pairs <- list()
   vlength <- 0
   iter <- 0
   remaining_idx <- seq_len(n_combs)
   
   while (vlength < n_pairs && length(remaining_idx) > 1 && iter < max_iter) {
      iter <- iter + 1
      message("  Iteration ", iter, ": sampling from ", length(remaining_idx), " remaining combinations.")
      
      # Sample a manageable chunk
      chunk_idx <- sample(remaining_idx, size = min(chunk_size, length(remaining_idx)))
      remaining_idx <- setdiff(remaining_idx, chunk_idx)
      
      sampled_combs <- combs_k[chunk_idx]
      sampled_bits  <- bitsets[chunk_idx]
      n_sampled <- length(sampled_combs)
      
      # Parallel disjoint search (lightweight)
      valid_chunk <- bplapply(seq_len(n_sampled - 1), function(i) {
         bit_i <- sampled_bits[[i]]
         # Compare against later combinations efficiently
         overlaps <- vapply((i + 1):n_sampled, function(j) any(bit_i & sampled_bits[[j]]), logical(1))
         disjoint_js <- which(!overlaps) + i
         if (length(disjoint_js) == 0) return(NULL)
         lapply(disjoint_js, function(j) list(A = dss[sampled_combs[[i]]], B = dss[sampled_combs[[j]]]))
      }, BPPARAM = bp_param)
      
      valid_chunk <- unlist(valid_chunk, recursive = FALSE)
      message("    Found ", length(valid_chunk), " disjoint pairs in this chunk.")
      
      if (length(valid_chunk) > 0) {
         valid_pairs <- c(valid_pairs, valid_chunk)
         vlength <- length(valid_pairs)
      }
      
      if(length(valid_chunk) < 500){
         chunk_size <- chunk_size * 2
      }
      
      if (vlength >= n_pairs) {
         message("  Reached target of ", n_pairs, " disjoint pairs.")
         valid_pairs <- valid_pairs[seq_len(n_pairs)]
         break
      }
   }
   
   if (iter >= max_iter)
      message("  Stopped early after reaching max iterations (", max_iter, ").")
   
   message("Total disjoint pairs collected for size ", k, ": ", vlength)
   
   if (vlength == 0L) {
      all_results[[k]] <- data.frame(GroupSize = k, Correlation = NA_real_)
      next
   }
   
   # Compute correlations in parallel
   message("Computing rankings for ", vlength, " disjoint pairs for size ", k)
   results_k <- bplapply(valid_pairs, BPPARAM = bp_param, function(pair) {
      rank_A <- compute_ranking(pair$A, alldf)
      rank_B <- compute_ranking(pair$B, alldf)
      
      common <- intersect(rank_A, rank_B)
      if (length(common) < 2) return(NULL)
      
      r1 <- match(common, rank_A)
      r2 <- match(common, rank_B)
      
      data.frame(
         Correlation = suppressWarnings(cor(r1, r2, method = "spearman")),
         GroupSize = k,
         Comparison = paste(paste(pair$A, collapse = "-"), paste(pair$B, collapse = "-"), sep = "_vs_")
      )
   })
   
   all_results[[k]] <- data.table::rbindlist(results_k, fill = TRUE)
}

final_df <- data.table::rbindlist(all_results, fill = TRUE)
output_dir <- file.path(savedir, paste0(type, ".rds"))
saveRDS(final_df, output_dir)

message_time("Completed pipeline")
