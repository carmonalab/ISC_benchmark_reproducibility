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
   
   cor_test <- cor.test(mean_metrics$mean_f1, mean_metrics$mean_isc)
   stats <- data.frame(
      dataset = dataset_label, isc = isc_col, isc_label = isc_label, group_by = "classifier",
      r2 = unname(cor_test$estimate^2), pval = cor_test$p.value, n = nrow(mean_metrics)
   )
   
   p <- ggplot(mean_metrics, aes(y = mean_f1, x = mean_isc)) +
      geom_smooth(method = "lm", color = "black", linetype = "dashed", se = TRUE, alpha = 0.2) +
      geom_point(aes(color = classifier), size = 2, show.legend = FALSE) +
      ggrepel::geom_label_repel(aes(label = classifier, color = classifier),
                                size = 2, show.legend = FALSE, max.overlaps = 20) +
      scale_color_manual(values = palette) +
      labs(
         title = gsub("_", "-", dataset_label),
         subtitle = sprintf("R² = %.3f, p-value = %.2e", stats$r2, stats$pval),
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
      summarise(
         r2 = cor(mean_f1, score, use = "complete.obs")^2,
         pval = cor.test(mean_f1, score)$p.value,
         .groups = "drop"
      ) %>%
      mutate(dataset = dataset_label)
   
   legend_labels <- setNames(
      sprintf("%s\nR²=%.2f\np=%.1e", corrs$metric, corrs$r2, corrs$pval),
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
plot_ranking_comparison <- function(df, isc_col, palette, excl = "ground_truth", dataset_label = "") {
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
      scale_color_manual(values = palette) +
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
plot_ranking_per_dataset <- function(df, isc_col, palette, excl = "ground_truth",
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
      scale_color_manual(values = palette) +
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
plot_isc_vs_f1_grouped <- function(df, isc_col, isc_label, group_by = c("classifier", "cell_type"),
                                   palette = NULL, excl = "ground_truth", dataset_label = "") {
   group_by <- match.arg(group_by)
   group_cols <- c("dataset_id", group_by)
   
   mean_metrics <- df %>%
      filter(!is.na(f1) & !is.na(.data[[isc_col]])) %>%
      filter(!classifier %in% excl) %>%
      group_by(across(all_of(group_cols))) %>%
      summarise(
         mean_f1 = mean(f1, na.rm = TRUE),
         mean_isc = mean(.data[[isc_col]], na.rm = TRUE),
         .groups = "drop"
      )
   
   cor_test <- cor.test(mean_metrics$mean_f1, mean_metrics$mean_isc)
   stats <- data.frame(
      dataset = dataset_label, isc = isc_col, isc_label = isc_label, group_by = group_by,
      r2 = unname(cor_test$estimate^2), pval = cor_test$p.value, n = nrow(mean_metrics)
   )
   
   p <- ggplot(mean_metrics, aes(y = mean_f1, x = mean_isc, color = .data[[group_by]])) +
      geom_point(alpha = 0.6, size = 1.5) +
      geom_smooth(method = "lm", color = "black", linetype = "dashed", se = TRUE) +
      scale_x_continuous(limits = c(-0.1, 1)) +
      scale_y_continuous(limits = c(-0.1, 1)) +
      labs(
         title = dataset_label,
         subtitle = sprintf("R² = %.3f, p-value = %.2e", stats$r2, stats$pval),
         y = "F1 Score Test",
         x = paste0("ISC Test (", isc_label, ")"),
         color = tools::toTitleCase(gsub("_", " ", group_by))
      ) +
      ggpubr::theme_classic2() +
      theme(legend.position = "right")
   
   if (!is.null(palette) && group_by == "classifier") p <- p + scale_color_manual(values = palette)
   
   list(plot = p, stats = stats)
}
