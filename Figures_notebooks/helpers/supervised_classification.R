collapse_replicates <- function(x) {
   x %>%
      group_by(dataset_id, classifier, cell_type) %>%
      summarise(
         across(
            c(f1, accuracy, all_of(iscs)),
            ~ mean(.x, na.rm = TRUE)
         ),
         n_replicates = n_distinct(replicate),
         .groups = "drop"
      )
}


# Generate ranking-based color palette
rankings_palette <- function(cons_combined,
                             isc = "product",
                             verbose = FALSE)
{
   rankings_for_palette <- cons_combined %>%
      filter(!is.na(f1) & !is.na(.data[[isc]]),
             classifier != "ground_truth") %>%
      group_by(dataset_id, classifier) %>%
      summarise(
         mean_product_per_dataset = mean(.data[[isc]], na.rm = TRUE),
         .groups = "drop"
      ) %>%
      group_by(classifier) %>%
      summarise(
         geom_mean_product = exp(mean(log(mean_product_per_dataset), na.rm = TRUE)),
         mean_product = mean(mean_product_per_dataset, na.rm = T),
         .groups = "drop"
      ) %>% 
      mutate(rank = rank(-geom_mean_product, ties.method = "average")) %>%
      arrange(rank)
   
   if(verbose) {print(rankings_for_palette)}
   
   # Create a unified, elegant, color-blind friendly palette
   # Integrating highly distinguishable top 3 best and top 3 worst colors
   n_classifiers <- nrow(rankings_for_palette)
   
   # Define the 6 distinct colors for top 3 best and top 3 worst
   colors_top3_best <- c("#0173B2", "#029E73", "#B19CD9")  # Deep blue, forest green, malva/violet
   colors_top3_worst <- c("#DE8F05", "#FFD700", "#882255")  # Orange, gold, purple-red
   
   if (n_classifiers <= 6) {
      # If 6 or fewer classifiers, use the 6 distinct colors
      gradient_colors <- c(colors_top3_best, colors_top3_worst)[1:n_classifiers]
   } else {  
      # For more than 6 classifiers, create a palette that:
      # - Uses top 3 best colors for the best performers
      # - Interpolates intermediate colors
      # - Uses top 3 worst colors for the worst performers
      
      # Generate intermediate colors using color interpolation
      # From best (blue) to worst (red)
      intermediate_palette <- colorRampPalette(c("#0173B2", "#E0E0E0", "#882255"))(n_classifiers - 6)
      
      # Combine: top 3 best + intermediate + top 3 worst
      gradient_colors <- c(
         colors_top3_best,
         intermediate_palette[1:(n_classifiers - 6)],
         colors_top3_worst
      )  
   }
   
   # Create named palette
   classifier_colors <- setNames(gradient_colors, rankings_for_palette$classifier)
   
   # Identify top 6 classifiers
   top_3_best <- rankings_for_palette$classifier[1:3]
   top_3_worst <- tail(rankings_for_palette$classifier, 3)
   top_6_classifiers <- c(top_3_best, top_3_worst)
   
   # For f1_vs_product_mean plot: highlight top 6, gray out the rest
   classifier_colors_highlight <- classifier_colors
   classifier_colors_highlight[!names(classifier_colors_highlight) %in% top_6_classifiers] <- "#808080"
   
   return(classifier_colors)
}



# Pearson R2 can be inflated by a few points sitting near the ISC ceiling (close to 1), which act
# as high-leverage points on the linear fit. Report Spearman (rank-based, robust to this) alongside
# Pearson, and flag points with Cook's distance > 4/n (rule of thumb) as influential.
# Also report the observed range/SD of the ISC metric itself: a metric compressed into a narrow
# range (e.g. cLISI often sits within 0.75-1) has much less room to vary than one spanning the full
# 0-1 range (e.g. ISC Local), so R2 values are not directly comparable across metrics with very
# different observed ranges (restriction-of-range problem) - a narrow-range metric needs a much
# tighter fit to reach the same R2, making any high R2 there less robust/more leverage-sensitive.
cor_stats <- function(isc, f1, min_n = 3) {
   keep <- stats::complete.cases(isc, f1)
   isc <- isc[keep]; f1 <- f1[keep]
   n <- length(isc)
   if (n < min_n || stats::sd(isc) == 0 || stats::sd(f1) == 0) {
      return(data.frame(
         r2 = NA_real_, pval = NA_real_,
         r2_spearman = NA_real_, pval_spearman = NA_real_,
         n_influential = NA_integer_, n = n,
         isc_min = suppressWarnings(min(isc)), isc_max = suppressWarnings(max(isc)),
         isc_range = NA_real_, isc_sd = stats::sd(isc)
      ))
   }
   pear <- stats::cor.test(isc, f1, method = "pearson")
   spear <- suppressWarnings(stats::cor.test(isc, f1, method = "spearman"))
   cooks <- stats::cooks.distance(stats::lm(f1 ~ isc))
   data.frame(
      r2 = unname(pear$estimate^2), pval = pear$p.value,
      r2_spearman = unname(spear$estimate^2), pval_spearman = spear$p.value,
      n_influential = sum(cooks > 4 / n), n = n,
      isc_min = min(isc), isc_max = max(isc), isc_range = max(isc) - min(isc), isc_sd = stats::sd(isc)
   )
}

# Mixed-effects model: f1 ~ ISC + classifier + (1 | dataset_id). Unlike cor_stats()/
# compute_r2_per_dataset_celltype(), this adjusts for classifier identity and avoids
# pseudo-replication from the repeated classifier/cell-type observations within each dataset
# (random intercept per dataset_id). Only base lme4 is available in this project (no lmerTest/
# performance), so the ISC fixed-effect p-value uses a normal (Wald z) approximation instead of a
# Satterthwaite-corrected t-test, and r2_marginal/r2_conditional are pseudo-R2 (squared correlation
# between predicted and observed F1, fixed-effects-only vs fixed+random), not the exact
# Nakagawa & Schielzeth decomposition.
# The ISC metric is z-scored (mean 0, SD 1) before fitting: this avoids convergence issues from
# very different ISC scales/ranges, and makes `estimate` comparable across metrics (change in F1
# per 1 SD increase in that ISC) rather than confounded by each metric's own observed range.
fit_isc_f1_lmer <- function(df, isc_col, isc_label = NULL, excl = "ground_truth",
                            model = ~ classifier + (1 | dataset_id) + (1 | dataset_id:replicate)) {
   if (is.null(isc_label)) isc_label <- labs[match(isc_col, iscs)]

   dat <- df %>%
      filter(!is.na(f1) & !is.na(.data[[isc_col]])) %>%
      filter(!classifier %in% excl) %>%
      mutate(classifier = droplevels(factor(classifier)))

   n_datasets <- dplyr::n_distinct(dat$dataset_id)
   n_classifiers <- dplyr::n_distinct(dat$classifier)
   isc_mean <- mean(dat[[isc_col]])
   isc_sd <- stats::sd(dat[[isc_col]])

   empty <- data.frame(
      isc = isc_col, isc_label = isc_label, estimate = NA_real_, se = NA_real_,
      ci_low = NA_real_, ci_high = NA_real_, pval = NA_real_,
      r2_marginal = NA_real_, r2_conditional = NA_real_,
      n = nrow(dat), n_datasets = n_datasets, n_classifiers = n_classifiers, singular = NA
   )
   if (nrow(dat) < 10 || n_datasets < 2 || n_classifiers < 2 || isc_sd == 0) return(empty)

   dat$.isc_z <- (dat[[isc_col]] - isc_mean) / isc_sd
   model_formula <- stats::update.formula(model, f1 ~ .isc_z + .)

   fit <- tryCatch(
      lme4::lmer(model_formula,
                 data = dat, REML = TRUE,
                control = lme4::lmerControl(check.conv.singular = "ignore")),
      error = function(e) NULL
   )
   if (is.null(fit) || !".isc_z" %in% rownames(summary(fit)$coefficients)) return(empty)

   co <- summary(fit)$coefficients
   est <- unname(co[".isc_z", "Estimate"])
   se  <- unname(co[".isc_z", "Std. Error"])
   pval <- 2 * stats::pnorm(-abs(est / se))

   r2_marginal <- suppressWarnings(stats::cor(stats::predict(fit, re.form = NA), dat$f1)^2)
   r2_conditional <- suppressWarnings(stats::cor(stats::predict(fit), dat$f1)^2)

   data.frame(
      isc = isc_col, isc_label = isc_label,
      estimate = est, se = se, ci_low = est - 1.96 * se, ci_high = est + 1.96 * se, pval = pval,
      r2_marginal = r2_marginal, r2_conditional = r2_conditional,
      n = nrow(dat), n_datasets = n_datasets, n_classifiers = n_classifiers,
      singular = lme4::isSingular(fit)
   )
}

# Compress [-0.1, 0.7] into 25% of the axis, [0.7, 1] into the remaining 75%
compressed_axis_trans <- scales::trans_new(
   name = "compressed",
   transform = function(x) {
      ifelse(x <= 0.8,
             (x + 0.1) / 3.2,           # [-0.1, 0.7] → [0, 0.25]
             0.25 + (x - 0.8) * 2.5)    # [0.7, 1] → [0.25, 1]
   },
   inverse = function(x) {
      ifelse(x <= 0.25,
             x * 3.2 - 0.1,             # [0, 0.25] → [-0.1, 0.7]
             (x - 0.25) / 2.5 + 0.8)    # [0.25, 1] → [0.7, 1]
   }
)

# Scatter of mean F1 vs mean ISC per classifier, labelled points
plot_isc_vs_f1_per_classifier <- function(df,
                                          isc_col,
                                          isc_label = NULL,
                                          palette = "f1",
                                          excl = "ground_truth",
                                          dataset_label = "") {
   palette <- rankings_palette(df, palette)
   if(is.null(isc_label)){
      isc_label <- labs[match(isc_col, iscs)]
   }
   
   mean_metrics <- df %>%
      filter(!is.na(f1) & !is.na(.data[[isc_col]])) %>%
      filter(!classifier %in% excl) %>%
      mutate(classifier = factor(classifier, levels = names(palette))) %>%
      group_by(classifier) %>%
      summarise(
         mean_f1 = mean(f1, na.rm = TRUE),
         mean_isc = mean(.data[[isc_col]], na.rm = TRUE),
         .groups = "drop"
      )
   
   stats <- cbind(
      data.frame(dataset = dataset_label, isc = isc_col, isc_label = isc_label, group_by = "classifier"),
      cor_stats(mean_metrics$mean_isc, mean_metrics$mean_f1)
   )
   
   p <- ggplot(mean_metrics, aes(y = mean_f1, x = mean_isc)) +
      geom_smooth(method = "lm", color = "black", linetype = "dashed", se = TRUE, alpha = 0.2) +
      geom_point(aes(color = classifier), size = 2, show.legend = FALSE) +
      ggrepel::geom_label_repel(aes(label = classifier, color = classifier),
                                size = 2, show.legend = FALSE, max.overlaps = 20) +
      scale_color_manual(values = palette) +
      labs(
         title = gsub("_", "-", dataset_label),
         subtitle = sprintf("R²(Pearson) = %.3f (p = %.2e) | R²(Spearman) = %.3f (p = %.2e)%s",
                            stats$r2, stats$pval, stats$r2_spearman, stats$pval_spearman,
                            ifelse(stats$n_influential > 0, sprintf(" | %d influential pt(s)", stats$n_influential), "")),
         y = "F1 Score Test",
         x = paste0("ISC Test (", isc_label, ")"),
         color = "Classifier"
      ) +
      ggpubr::theme_classic2()
   
   list(plot = p, stats = stats)
}

# Overlay of all ISC metrics vs F1 (per-dataset mean, one fit line per metric)
plot_multi_isc_corr <- function(df, iscs, labs, isc_colors, excl = "ground_truth", dataset_label = "") {
   mean_metrics <- df %>%
      filter(!classifier %in% excl) %>%
      group_by(dataset_id, classifier) %>%
      summarise(
         mean_f1 = mean(f1, na.rm = TRUE),
         across(all_of(iscs), ~ mean(.x, na.rm = TRUE), .names = "{.col}"),
         .groups = "drop"
      ) %>%
      pivot_longer(-c(dataset_id, classifier, mean_f1), names_to = "metric", values_to = "score") %>%
      mutate(metric = factor(metric, levels = iscs, labels = labs))
   
   corrs <- mean_metrics %>%
      group_by(metric) %>%
      group_modify(~ cor_stats(.x$score, .x$mean_f1)) %>%
      ungroup() %>%
      mutate(dataset = dataset_label)
   
   # Pearson R2 next to Spearman R2 so ceiling-effect inflation (e.g. ISC metrics clustered near 1)
   # is visible directly in the legend
   legend_labels <- setNames(
      sprintf("%s\nR²p=%.2f R²s=%.2f\np=%.1e", corrs$metric, corrs$r2, corrs$r2_spearman, corrs$pval),
      corrs$metric
   )
   
   p <- ggplot(mean_metrics, aes(y = mean_f1, x = score, color = metric, fill = metric)) +
      geom_smooth(method = "lm", se = TRUE, linetype = "solid", alpha = 0.4, show.legend = TRUE) +
      scale_color_manual(values = isc_colors, labels = legend_labels, name = "ISC metric") +
      scale_fill_manual(values = isc_colors, guide = "none") +
      labs(title = dataset_label, y = "F1 Score Test", x = "ISC Test") +
      guides(fill = "none") +
      ggpubr::theme_classic2() +
      theme(legend.position = "bottom", legend.title = element_blank())
   
   list(plot = p, stats = corrs)
}

# Overall ranking agreement between F1 and a given ISC metric (geometric mean across datasets)
plot_ranking_comparison <- function(df, isc_col, excl = "ground_truth", dataset_label = "") {
   rankings <- df %>%
      filter(!is.na(f1) & !is.na(.data[[isc_col]])) %>%
      filter(!classifier %in% excl) %>%
      group_by(dataset_id, classifier) %>%
      summarise(
         mean_f1_per_dataset = mean(f1, na.rm = TRUE),
         mean_isc_per_dataset = mean(.data[[isc_col]], na.rm = TRUE),
         .groups = "drop"
      ) %>%
      group_by(classifier) %>%
      summarise(
         geom_mean_f1 = exp(mean(log(mean_f1_per_dataset), na.rm = TRUE)),
         geom_mean_isc = exp(mean(log(mean_isc_per_dataset), na.rm = TRUE)),
         .groups = "drop"
      ) %>%
      mutate(
         rank_f1 = rank(-geom_mean_f1, ties.method = "average"),
         rank_isc = rank(-geom_mean_isc, ties.method = "average"),
         rank_difference = abs(rank_f1 - rank_isc)
      ) %>%
      mutate(classifier = factor(classifier, levels = names(palette))) %>%
      arrange(rank_isc)
   
   rank_cor <- cor.test(rankings$rank_f1, rankings$rank_isc, method = "spearman")
   
   p <- ggplot(rankings, aes(x = rank_isc, y = rank_f1)) +
      geom_point(aes(color = classifier), size = 4) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40") +
      ggrepel::geom_label_repel(aes(label = classifier, color = classifier), size = 5) +
      scale_x_reverse(breaks = seq_len(nrow(rankings))) +
      scale_y_reverse(breaks = seq_len(nrow(rankings))) +
      labs(
         title = "Ranking Agreement: ISC vs F1",
         subtitle = dataset_label,
         x = "ISC Rank", y = "F1 Rank", color = "Classifier"
      ) +
      coord_fixed() +
      ggpubr::theme_classic2() +
      theme(legend.position = "none")
   
   list(plot = p, rankings = rankings, rank_cor = rank_cor)
}

# Per-dataset ranking of classifiers (ranked within each dataset, then averaged across datasets)
plot_ranking_per_dataset <- function(df, isc_col, excl = "ground_truth",
                                     use_geom_mean = TRUE, dataset_label = "") {
   dataset_metrics_mean <- df %>%
      filter(!is.na(f1) & !is.na(.data[[isc_col]])) %>%
      filter(!classifier %in% excl) %>%
      group_by(dataset_id, classifier) %>%
      summarise(
         mean_f1 = mean(f1, na.rm = TRUE),
         mean_isc = mean(.data[[isc_col]], na.rm = TRUE),
         .groups = "drop"
      )
   
   dataset_rankings <- dataset_metrics_mean %>%
      group_by(dataset_id) %>%
      mutate(
         rank_f1 = rank(-mean_f1, ties.method = "average"),
         rank_isc = rank(-mean_isc, ties.method = "average")
      ) %>%
      ungroup()
   
   mean_rank_per_classifier <- dataset_rankings %>%
      group_by(classifier) %>%
      summarise(
         mean_rank_f1 = if (use_geom_mean) exp(mean(log(rank_f1), na.rm = TRUE)) else mean(rank_f1, na.rm = TRUE),
         mean_rank_isc = if (use_geom_mean) exp(mean(log(rank_isc), na.rm = TRUE)) else mean(rank_isc, na.rm = TRUE),
         sd_rank_f1 = sd(rank_f1, na.rm = TRUE),
         sd_rank_isc = sd(rank_isc, na.rm = TRUE),
         n_datasets = n(),
         .groups = "drop"
      ) %>%
      mutate(classifier = factor(classifier, levels = names(palette))) %>%
      arrange(mean_rank_isc)
   
   rank_cor_dataset <- cor.test(mean_rank_per_classifier$mean_rank_f1,
                                mean_rank_per_classifier$mean_rank_isc, method = "spearman")
   mean_type <- if (use_geom_mean) "Geometric Mean" else "Arithmetic Mean"
   
   p <- ggplot(mean_rank_per_classifier, aes(x = mean_rank_isc, y = mean_rank_f1)) +
      geom_point(aes(color = classifier), size = 4) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40") +
      ggrepel::geom_label_repel(aes(label = classifier, color = classifier), size = 5) +
      scale_x_reverse() +
      scale_y_reverse() +
      labs(
         title = dataset_label,
         subtitle = sprintf("Ranking Agreement per Dataset: ISC vs F1 [%s], Spearman ρ = %.3f (p = %.2e)",
                            mean_type, rank_cor_dataset$estimate, rank_cor_dataset$p.value),
         x = "Mean ISC Rank", y = "Mean F1 Rank", color = "Classifier"
      ) +
      coord_fixed() +
      ggpubr::theme_classic2() +
      theme(legend.position = "none")
   
   list(plot = p, mean_rank_per_classifier = mean_rank_per_classifier, rank_cor = rank_cor_dataset)
}

# Generic scatter of mean F1 vs mean ISC, grouped/colored by either "classifier" or "cell_type"
plot_isc_vs_f1_grouped <- function(df,
                                   isc_col,
                                   isc_label,
                                   group_by = c("classifier", "cell_type"),
                                   palette = "f1",
                                   excl = "ground_truth",
                                   dataset_label = "") {
   
   group_cols <- c("dataset_id", group_by)
   
   palette <- rankings_palette(df, palette)
   
   mean_metrics <- df %>%
      filter(!is.na(f1) & !is.na(.data[[isc_col]])) %>%
      filter(!classifier %in% excl) %>%
      group_by(across(all_of(group_cols))) %>%
      summarise(
         mean_f1 = mean(f1, na.rm = TRUE),
         mean_isc = mean(.data[[isc_col]], na.rm = TRUE),
         .groups = "drop"
      )
   
   stats <- cbind(
      data.frame(
         dataset = dataset_label,
         isc = isc_col,
         isc_label = isc_label,
         group_by = paste(group_cols, collapse = ":")
      ),
      cor_stats(mean_metrics$mean_isc, mean_metrics$mean_f1)
   )
   
   p <- ggplot(mean_metrics, aes(y = mean_f1, x = mean_isc, color = classifier)) +
      geom_point(alpha = 0.6, size = 1.5) +
      geom_smooth(method = "lm", color = "black", linetype = "dashed", se = TRUE) +
      scale_x_continuous(limits = c(-0.1, 1)) +
      scale_y_continuous(limits = c(-0.1, 1)) +
      scale_color_manual(values = palette) +
      labs(
         title = dataset_label,
         subtitle = sprintf("R²(Pearson) = %.3f | R²(Spearman) = %.3f (p = %.2e)%s",
                            stats$r2, stats$r2_spearman, stats$pval_spearman,
                            ifelse(stats$n_influential > 0, sprintf(" | %d influential pt(s)", stats$n_influential), "")),
         y = "F1 Score Test",
         x = paste0("ISC Test (", isc_label, ")"),
         color = tools::toTitleCase(gsub("_", " ", group_by))
      ) +
      ggpubr::theme_classic2() +
      theme(legend.position = "right")

   
   list(plot = p, stats = stats)
}

# Per-dataset R2 (F1 vs ISC) computed at the cell-type level (one fit per dataset_id x replicate,
# using classifier x cell_type pairs as data points), so variability across replicates is preserved
compute_r2_per_dataset_celltype <- function(df,
                                            isc_col,
                                            group_cols = c("dataset_id"),
                                            excl = "ground_truth",
                                            min_n = 3) {
   df %>%
      filter(!is.na(f1) & !is.na(.data[[isc_col]])) %>%
      filter(!classifier %in% excl) %>%
      group_by(across(all_of(group_cols))) %>%
      filter(n() >= min_n) %>%
      group_modify(~ cor_stats(.x[[isc_col]], .x$f1, min_n = min_n)) %>%
      ungroup() %>%
      mutate(isc = isc_col)
}

# Join prediction F1 (non-ground_truth, non-excluded classifiers) with each ISC metric computed on
# the ground_truth reference row, matched by dataset_id/cell_type. Loops over every metric in
# isc_col (defaults to the global `iscs`, skipping classifiers in excl) and returns a named list of
# joined data frames, each keeping its own isc_col column name so it can be fed directly into
# cor_stats()/compute_r2_per_dataset_celltype()/fit_isc_f1_lmer().
join_reference_consistency <- function(df, isc_col = iscs, excl = "ground_truth") {
   pred_summary <- df %>%
      filter(classifier != "ground_truth", !classifier %in% excl, !is.na(f1)) %>%
      select(cell_type, dataset_id, f1, classifier)

   lapply(isc_col, function(ic) {
      cons_ref <- df %>%
         filter(classifier == "ground_truth", !is.na(.data[[ic]])) %>% 
         select(-f1, -classifier)

      pred_summary %>%
         left_join(cons_ref, by = c("dataset_id", "cell_type")) %>%
         filter(!is.na(.data[[ic]]))
   }) %>%
      setNames(isc_col)
}

# Quadrant contingency of reference ISC vs prediction F1, thresholded at ths_isc/ths_f1: computes
# per-quadrant counts/percentages and renders the bin2d density plot with quadrant labels. `joined`
# is expected to come from join_reference_consistency()[[isc_col]].
plot_isc_f1_contingency <- function(joined, isc_col, isc_label = NULL,
                                    ths_isc = 0.5, ths_f1 = 0.5, dataset_label = "") {
   if (is.null(isc_label)) isc_label <- labs[match(isc_col, iscs)]

   quadrant_counts <- joined %>%
      mutate(quadrant = case_when(
         .data[[isc_col]] >= ths_isc & f1 >= ths_f1 ~ "top-right",
         .data[[isc_col]] <  ths_isc & f1 >= ths_f1 ~ "top-left",
         .data[[isc_col]] >= ths_isc & f1 <  ths_f1 ~ "bottom-right",
         TRUE                                       ~ "bottom-left"
      )) %>%
      count(quadrant) %>%
      mutate(
         dataset = dataset_label, isc = isc_col, isc_label = isc_label,
         prop = n / sum(n),
         x = case_when(quadrant == "top-right"   ~ 0.75,
                       quadrant == "top-left"     ~ 0.25,
                       quadrant == "bottom-right" ~ 0.75,
                       quadrant == "bottom-left"  ~ 0.25),
         y = case_when(quadrant == "top-right"   ~ 0.85,
                       quadrant == "top-left"     ~ 0.85,
                       quadrant == "bottom-right" ~ 0.25,
                       quadrant == "bottom-left"  ~ 0.25),
         label = paste0(sprintf("%.1f", 100 * prop), "%")
      )

   p <- joined %>%
      ggplot(aes(x = .data[[isc_col]], y = f1)) +
      geom_bin2d(bins = 50, show.legend = TRUE) +
      scale_fill_viridis_c(option = "magma", trans = "sqrt") +
      geom_vline(xintercept = ths_isc, linetype = "dashed", color = "#00FFFF", linewidth = 0.8) +
      geom_hline(yintercept = ths_f1, linetype = "dashed", color = "#00FFFF", linewidth = 0.8) +
      geom_label(data = quadrant_counts,
                aes(x = x, y = y, label = label),
                inherit.aes = FALSE, size = 6, fill = "grey90", color = "black", alpha = 0.9) +
      # scale_x_continuous(limits = c(0, 1)) +
      # scale_y_continuous(limits = c(0, 1)) +
      labs(
         title = paste0(gsub("_", "-", dataset_label), " - ", isc_label),
         subtitle = "Reference ISC vs Prediction F1, per cell type",
         x = paste0("Reference ISC (", isc_label, ")"),
         y = "F1 score Test"
      ) +
      ggpubr::theme_classic2() +
      theme(legend.position = "right")

   list(plot = p, quadrant_counts = quadrant_counts)
}


plot_iscTest_vs_f1 <- function(df,
                               isc_col,
                               isc_label,
                               excl = "ground_truth",
                               dataset_label = "") {

   stats <- cbind(
      data.frame(
         dataset = dataset_label,
         isc = isc_col,
         isc_label = isc_label
      ),
      cor_stats(df$f1, df[[isc_col]])
   )
   
   p <- ggplot(mean_metrics, aes(y = mean_f1, x = mean_isc, color = classifier)) +
      geom_point(alpha = 0.6, size = 1.5) +
      geom_smooth(method = "lm", color = "black", linetype = "dashed", se = TRUE) +
      scale_x_continuous(limits = c(-0.1, 1)) +
      scale_y_continuous(limits = c(-0.1, 1)) +
      scale_color_manual(values = palette) +
      labs(
         title = dataset_label,
         subtitle = sprintf("R²(Pearson) = %.3f | R²(Spearman) = %.3f (p = %.2e)%s",
                            stats$r2, stats$r2_spearman, stats$pval_spearman,
                            ifelse(stats$n_influential > 0, sprintf(" | %d influential pt(s)", stats$n_influential), "")),
         y = "F1 Score Test",
         x = paste0("ISC Test (", isc_label, ")"),
         color = tools::toTitleCase(gsub("_", " ", group_by))
      ) +
      ggpubr::theme_classic2() +
      theme(legend.position = "right")
   
   
   list(plot = p, stats = stats)
}
