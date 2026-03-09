suppressPackageStartupMessages({
  library(dplyr)
})

optimal_dataset <- function(obj,
                            ident = "celltype",
                            sample = "sample",
                            min_samples = 9,
                            max_samples = 15,
                            min_cells_sample = 200,
                            verbose = TRUE) {
  md <- obj@meta.data %>%
    mutate(ident = .data[[ident]],
           sample = .data[[sample]]) %>%
    select(sample, ident)

  nsamples <- length(unique(md[["sample"]]))

  nidents <- md %>%
    group_by(sample, ident) %>%
    summarize(ncells_ident = n(), .groups = "drop") %>%
    group_by(sample) %>%
    mutate(nident = n(),
           ncells_sample = sum(ncells_ident))

  if (nsamples > max_samples) {
    top_rep <- nidents %>%
      distinct(sample, nident, ncells_sample) %>%
      filter(nident > 1 & ncells_sample >= min_cells_sample) %>%
      arrange(desc(nident), desc(ncells_sample)) %>%
      head(max_samples) %>%
      pull(sample)

    obj <- obj[, obj@meta.data[[sample]] %in% top_rep]
    if (verbose) {
      message("- Subsetting dataset from ", nsamples, " to ", length(top_rep), " samples")
    }
    nsamples <- length(unique(obj@meta.data[[sample]]))
  }

  if (nsamples < min_samples) {
    if (verbose) message("- Only ", nsamples, " samples; returning NULL")
    return(NULL)
  }

  if (verbose) message("- Returning ", nsamples, " samples")
  obj
}

split_samples_query_reference <- function(metadata,
                                         sample_col = "sample",
                                         prop_query = 0.7,
                                         seed = 22) {
  allsamples <- unique(metadata[[sample_col]])
  nsamples <- floor(length(allsamples) * prop_query)
  set.seed(seed)
  query_samples <- sample(allsamples, size = nsamples)
  list(query_samples = query_samples)
}
