library(dplyr)

perturb_cells <- function(counts,
                          cells_perturb,
                          cells_reference,
                          noise = 0.5,
                          seed = 22) {
  stopifnot(all(c(cells_perturb, cells_reference) %in% colnames(counts)))
  stopifnot(is.numeric(noise) && length(noise) == 1 && noise >= 0)
  counts <- as.matrix(counts)
  set.seed(seed) # for reproducibility
  # Vectorized perturbation using apply
  counts[, cells_perturb] <- apply(counts[, cells_perturb, drop = FALSE], 2, function(original) {
    # take one random cell from reference
    ref <- sample(cells_reference, 1)
    # difference of expression
    diff <- original - counts[, ref]
    # apply rate
    new_counts <- original + diff * noise
    new_counts[new_counts < 0] <- 0
    return(new_counts)
  })
  return(counts)
}


get_consistency_wide <- function(scTypeEval,
                                 consistency_metric = c("silhouette", "2label_silhouette"),
                                 global = "2label_silhouette | Pseudobulk:Cosine",
                                 local = "silhouette | recip_classif:Match",
                                 wide = TRUE){
   
   consis <- get_consistency(scTypeEval,
                             consistency_metric = consistency_metric)
   consis <- consis %>% 
      mutate(method_type = paste(consistency_metric,
                                 dissimilarity_method,
                                 sep = " | ")) %>% 
      filter(method_type %in% c(global, local)) %>% 
      mutate(method_type = factor(method_type,
                                  levels = c(global, local),
                                  labels = c("global", "local")))
   if(wide){
   consis <- consis %>% 
      select(method_type, measure, celltype) %>% 
      pivot_wider(names_from = "method_type",
                  values_from = "measure") %>% 
      mutate(product = global * local)
   }
   return(consis)
}


select_perturbed_cells <- function(metadata,
                                   ident = "OriginalAnnotationLevel1",
                                   sample = "sample",
                                   perturb_celltype = "B_cell",
                                   type = c("per_sample", # perturb whole different samples
                                              "same_sample", # perturb same ratio in each sample
                                              "different_sample" # perturb different ratio per sample
                                            ),
                                   rate_same_sample = 0.5,
                                   rates_different_samples = 0.25,
                                   seed = 22){
  # get only needed cols
  metadata$ident <- metadata[, ident]
  metadata$sample <- metadata[, sample]
  md <- metadata %>% 
    select(ident, sample)
  md$pertub_annot <- as.character(md$ident)
  
  md_perturb <- md %>% 
    filter(ident == perturb_celltype)
  
  unique_samples <- unique(metadata$sample)
  n_samples <- length(unique_samples)
  
  set_prop <- function(num = 0, rates_different_samples){
    if(num %% 2 == 0){
      rates_different_samples
    } else {
      1 - rates_different_samples
    }
  }
  
  set.seed(seed)
  cells_perturb <- switch (type,
                           "per_sample" = {
                             samples_perturb <- sample(unique_samples, floor(n_samples/2))
                             md_perturb %>% 
                               filter(sample %in% samples_perturb) %>% 
                               rownames()
                           },
                           "same_sample" = {
                             unlist(
                               lapply(split(rownames(md_perturb), md_perturb$sample), function(x) {
                                 sample(x, ceiling(length(x) * rate_same_sample))
                               })
                             )
                           },
                           "different_sample" = {
                             lapply(unique_samples, function(s){
                               cellsids <- md_perturb %>% 
                                 filter(sample == s) %>% 
                                 rownames()
                               pos <- which(unique_samples == s)
                               prop <- set_prop(pos, rates_different_samples)
                               ncells <- ceiling(prop * length(cellsids))
                               sample(cellsids, ncells)
                             }) %>% 
                               unlist()
                           }
  )
  
  # add metadata
  md[cells_perturb, "pertub_annot"] <- "perturbed"
  
  return(md)
}
