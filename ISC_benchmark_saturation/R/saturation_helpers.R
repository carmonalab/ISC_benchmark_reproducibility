#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(stringr)
  library(data.table)
  library(BiocParallel)
})

required_cols <- c(
  "dissimilarity_method",
  "consistency.metric",
  "method_type",
  "assay",
  "dataset",
  "ident",
  "score"
)

compute_ranking <- function(selected_dss, alldf, weight_vector) {
  alldf %>%
    mutate(dss = paste(dataset, ident, sep = "_")) %>%
    filter(dss %in% selected_dss) %>%
    group_by(method_type, assay) %>%
    summarise(score = mean(score, na.rm = TRUE), .groups = "drop") %>%
    group_by(method_type) %>%
    summarise(score = weighted.mean(score, w = unname(weight_vector[assay]), na.rm = TRUE), .groups = "drop") %>%
    mutate(score = round(score, 3)) %>%
    arrange(desc(score)) %>%
    pull(method_type)
}

load_merged_isc_input <- function(cfg) {
  merged_file <- cfg$paths$merged_results_file
  if (file.exists(merged_file)) {
    readRDS(merged_file)
  } else {
    NULL
  }
}

collect_metrics_from_dir <- function(metrics_dir) {
  if (!dir.exists(metrics_dir)) {
    stop("Fallback metrics directory does not exist: ", metrics_dir)
  }

  files <- list.files(
    metrics_dir,
    pattern = "_metrics\\.rds$",
    recursive = TRUE,
    full.names = TRUE
  )

  if (length(files) == 0) {
    stop("No *_metrics.rds files found under fallback directory: ", metrics_dir)
  }

  message(sprintf("Found %d metrics files for concatenation", length(files)))

  out <- lapply(files, function(f) {
    df <- readRDS(f)
    if (!is.data.frame(df)) {
      return(NULL)
    }
    df$source_file <- f
    df
  })

  out <- out[!vapply(out, is.null, logical(1))]
  if (length(out) == 0) {
    stop("No readable data frames in *_metrics.rds files")
  }

  data.table::rbindlist(out, fill = TRUE)
}

normalize_saturation_input <- function(df, cfg) {
  if (!is.data.frame(df)) {
    stop("Saturation input must be a data frame")
  }

  rename_map <- c(
    "dissimilarity.method" = "dissimilarity_method",
    "consistency_metric" = "consistency.metric",
    "dataset_id" = "dataset"
  )

  for (old_name in names(rename_map)) {
    new_name <- rename_map[[old_name]]
    if (old_name %in% names(df) && !(new_name %in% names(df))) {
      names(df)[names(df) == old_name] <- new_name
    }
  }

  if (!("score" %in% names(df)) && ("degradation_score" %in% names(df))) {
    df$score <- df$degradation_score
  }

  if (!("assay" %in% names(df)) && ("task" %in% names(df))) {
    df$assay <- unname(cfg$assay_map[as.character(df$task)])
  }

  if (!("method_type" %in% names(df)) &&
      all(c("consistency.metric", "dissimilarity_method") %in% names(df))) {
    df$method_type <- paste(df$consistency.metric, "|", df$dissimilarity_method)
  }

  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0) {
    stop(
      "Input data is missing required columns for saturation: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  df <- df %>%
    mutate(
      score = as.numeric(score),
      dataset = as.character(dataset),
      ident = as.character(ident),
      assay = as.character(assay),
      method_type = as.character(method_type),
      dissimilarity_method = as.character(dissimilarity_method),
      consistency.metric = as.character(consistency.metric)
    ) %>%
    filter(!is.na(score), !is.na(dataset), !is.na(ident), !is.na(assay)) %>%
    select(all_of(required_cols))

  if (nrow(df) == 0) {
    stop("No valid rows remain after normalization")
  }

  df
}

prepare_saturation_input <- function(merged_input, cfg) {
  if (is.data.frame(merged_input)) {
    candidate <- merged_input
  } else if (is.list(merged_input) && length(merged_input) > 0 && all(vapply(merged_input, is.data.frame, logical(1)))) {
    candidate <- data.table::rbindlist(merged_input, fill = TRUE)
  } else {
    candidate <- NULL
  }

  if (is.null(candidate)) {
    candidate <- collect_metrics_from_dir(cfg$paths$fallback_metrics_dir)
  }

  normalize_saturation_input(candidate, cfg)
}

as_bit <- function(indices, n_total) {
  bit <- logical(n_total)
  bit[indices] <- TRUE
  bit
}

find_disjoint_pairs <- function(combs_k, dss, n_total, n_pairs, chunk_size, max_iter, bp_param) {
  n_combs <- length(combs_k)
  if (n_combs < 2) {
    return(list())
  }

  if (n_combs > 500) {
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

    chunk_idx <- sample(remaining_idx, size = min(chunk_size, length(remaining_idx)))
    remaining_idx <- setdiff(remaining_idx, chunk_idx)

    sampled_combs <- combs_k[chunk_idx]
    sampled_bits <- bitsets[chunk_idx]
    n_sampled <- length(sampled_combs)

    if (n_sampled < 2) {
      next
    }

    valid_chunk <- bplapply(seq_len(n_sampled - 1), function(i) {
      bit_i <- sampled_bits[[i]]
      overlaps <- vapply((i + 1):n_sampled, function(j) any(bit_i & sampled_bits[[j]]), logical(1))
      disjoint_js <- which(!overlaps) + i
      if (length(disjoint_js) == 0) {
        return(NULL)
      }
      lapply(disjoint_js, function(j) list(A = dss[sampled_combs[[i]]], B = dss[sampled_combs[[j]]]))
    }, BPPARAM = bp_param)

    valid_chunk <- unlist(valid_chunk, recursive = FALSE)

    if (length(valid_chunk) > 0) {
      valid_pairs <- c(valid_pairs, valid_chunk)
      vlength <- length(valid_pairs)
    }

    if (length(valid_chunk) < 500) {
      chunk_size <- chunk_size * 2
    }

    if (vlength >= n_pairs) {
      valid_pairs <- valid_pairs[seq_len(n_pairs)]
      break
    }
  }

  valid_pairs
}

comparison_id_from_pair <- function(pair) {
  side_a <- paste(sort(pair$A), collapse = "-")
  side_b <- paste(sort(pair$B), collapse = "-")
  sides <- sort(c(side_a, side_b))
  paste(sides[1], sides[2], sep = "_vs_")
}

get_group_cache_file <- function(group_size, cfg) {
  cache_dir <- if (!is.null(cfg$paths$cache_dir)) cfg$paths$cache_dir else file.path(cfg$paths$output_dir, "cache_groups")
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  file.path(cache_dir, sprintf("group_size_%02d.rds", as.integer(group_size)))
}

load_group_cache <- function(group_size, cfg) {
  cache_file <- get_group_cache_file(group_size, cfg)
  if (!file.exists(cache_file)) {
    return(data.frame())
  }

  cached <- readRDS(cache_file)
  if (!is.data.frame(cached) || nrow(cached) == 0) {
    return(data.frame())
  }

  if (!("ComparisonID" %in% names(cached)) && ("Comparison" %in% names(cached))) {
    cached$ComparisonID <- as.character(cached$Comparison)
  }

  cached
}

save_group_cache <- function(group_size, results_df, cfg) {
  cache_file <- get_group_cache_file(group_size, cfg)
  saveRDS(results_df, cache_file)
  invisible(cache_file)
}

run_saturation_group_size <- function(group_size, alldf, cfg) {
  dss <- unique(paste(alldf$dataset, alldf$ident, sep = "_"))
  n_total <- length(dss)

  if ((group_size * 2) > n_total) {
    return(data.frame(GroupSize = group_size, Correlation = NA_real_))
  }

  set.seed(cfg$seed + as.integer(group_size))

  cached_results <- load_group_cache(group_size, cfg)
  if (nrow(cached_results) > 0) {
    cached_results <- cached_results %>%
      filter(!is.na(Correlation), !is.na(ComparisonID)) %>%
      distinct(ComparisonID, .keep_all = TRUE)
  }

  n_cached <- nrow(cached_results)
  n_needed <- max(0L, as.integer(cfg$analysis$n_pairs) - n_cached)

  if (n_needed == 0L) {
    message(sprintf("Group size %s already cached (%d/%d). Skipping.",
                    group_size, n_cached, cfg$analysis$n_pairs))
    return(cached_results)
  }

  combs_k <- combn(seq_len(n_total), m = group_size, simplify = FALSE)
  bp_param <- BiocParallel::MulticoreParam(
    workers = cfg$analysis$n_cores,
    progressbar = FALSE
  )

  valid_pairs <- find_disjoint_pairs(
    combs_k = combs_k,
    dss = dss,
    n_total = n_total,
    # Request enough pairs to go beyond the already-cached prefix,
    # otherwise deterministic sampling can repeatedly return only cached pairs.
    n_pairs = n_cached + max(n_needed * 3L, n_needed + 100L),
    chunk_size = cfg$analysis$chunk_size,
    max_iter = cfg$analysis$max_iter,
    bp_param = bp_param
  )

  if (length(valid_pairs) == 0) {
    if (nrow(cached_results) > 0) {
      return(cached_results)
    }
    return(data.frame(GroupSize = group_size, Correlation = NA_real_))
  }

  existing_ids <- if (nrow(cached_results) > 0) cached_results$ComparisonID else character(0)
  valid_pairs <- valid_pairs[!vapply(valid_pairs, function(p) comparison_id_from_pair(p) %in% existing_ids, logical(1))]

  if (length(valid_pairs) == 0) {
    return(cached_results)
  }

  if (length(valid_pairs) > n_needed) {
    valid_pairs <- valid_pairs[seq_len(n_needed)]
  }

  weight_vector <- unlist(cfg$weights)

  results_k <- bplapply(valid_pairs, BPPARAM = bp_param, function(pair) {
    comp_id <- comparison_id_from_pair(pair)

    rank_a <- compute_ranking(pair$A, alldf, weight_vector)
    rank_b <- compute_ranking(pair$B, alldf, weight_vector)

    common <- intersect(rank_a, rank_b)
    if (length(common) < 2) {
      return(NULL)
    }

    r1 <- match(common, rank_a)
    r2 <- match(common, rank_b)

    data.frame(
      Correlation = suppressWarnings(cor(r1, r2, method = "spearman")),
      GroupSize = group_size,
      ComparisonID = comp_id,
      Comparison = paste(paste(pair$A, collapse = "-"), paste(pair$B, collapse = "-"), sep = "_vs_"),
      stringsAsFactors = FALSE
    )
  })

  computed_results <- data.table::rbindlist(results_k, fill = TRUE)
  all_results <- data.table::rbindlist(list(cached_results, computed_results), fill = TRUE)

  if (nrow(all_results) > 0 && "ComparisonID" %in% names(all_results)) {
    all_results <- all_results %>%
      filter(!is.na(ComparisonID), !is.na(Correlation)) %>%
      distinct(ComparisonID, .keep_all = TRUE)
  }

  save_group_cache(group_size, all_results, cfg)
  all_results
}

save_saturation_outputs <- function(saturation_input, saturation_results, cfg) {
  out_dir <- cfg$paths$output_dir
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  input_file <- file.path(out_dir, sprintf("%s_input_concatenated.rds", cfg$analysis$type))
  results_file <- file.path(out_dir, sprintf("%s_saturation_results.rds", cfg$analysis$type))

  saveRDS(saturation_input, input_file)
  saturation_results_clean <- saturation_results
  if (is.data.frame(saturation_results_clean) && "ComparisonID" %in% names(saturation_results_clean)) {
    saturation_results_clean <- saturation_results_clean %>% select(-ComparisonID)
  }

  saveRDS(saturation_results_clean, results_file)

  list(
    input_file = input_file,
    results_file = results_file,
    n_input_rows = nrow(saturation_input),
    n_result_rows = nrow(saturation_results_clean)
  )
}
