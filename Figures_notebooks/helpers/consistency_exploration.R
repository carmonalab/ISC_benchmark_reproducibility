consistency_scatter <- function(df,
                                xval = "local",
                                yval = "global",
                                fill = "measure",
                                label = "celltype",
                                x_threshold = 0.45,
                                y_threshold = 0.74,
                                linewidth = 0.5,
                                alpha_line = 0.6,
                                xlabel = NULL,
                                ylabel = NULL) {
   
   if(is.null(xlabel)){
      xlabel <- paste0("ISC score (", xval, ")")
   }
   if(is.null(ylabel)){
      ylabel <- paste0("ISC score (", yval, ")")
   }
   
   
   pl_scatter <- df %>%
      ggplot(aes(x = .data[[xval]],
                 y = .data[[yval]])) +
      # threshold line of figure 5
      geom_vline(
         xintercept = x_threshold,
         linetype = "dashed",
         linewidth = linewidth,
         alpha = alpha_line
      ) +
      geom_hline(
         yintercept = y_threshold,
         linetype = "dashed",
         linewidth = linewidth,
         alpha = alpha_line
      ) +
      # Points
      geom_point(aes(fill = .data[[fill]]),
                 shape = 21,
                 size = 3.2,
                 stroke = 0.4,
                 color = "black",
                 alpha = 0.9,
                 show.legend = FALSE) +
      
      # Clean blue gradient (Nature-friendly, perceptually smooth)
      scale_fill_gradient(
         low  = "#F7FBFF",
         high = "#08306B"
      ) +
      
      # Repelled labels (better spacing control)
      ggrepel::geom_text_repel(
         aes(label = .data[[label]]),
         size = 3.5,
         max.overlaps = 8,
         box.padding = 0.4,
         point.padding = 0.3,
         segment.size = 0.3,
         segment.color = "grey50",
         min.segment.length = 0,
         seed = 123
      ) +
      
      labs(
         x = xlabel,
         y = ylabel
      ) +
      
      coord_cartesian(clip = "off") +
      
      theme_classic(base_size = 12) +
      theme(
         axis.title = element_text(size = 20),
         axis.text  = element_text(size = 11),
         axis.line  = element_line(linewidth = 0.4),
         axis.ticks = element_blank(),
         plot.margin = margin(10, 20, 10, 10)
      )
   
   return(pl_scatter)
}


consistency_barplot <- function(df,
                                xval = "local",
                                ylabel = "Author's annotation",
                                xlabel = NULL,
                                order = "celltype",
                                x_threshold = 0,
                                linewidth = 0.5,
                                alpha_line = 0.6) {
   
   if(is.null(xlabel)){
      xlabel <- paste0("ISC score (", xval, ")")
   }
   
   pl <- df %>% 
      ggplot(aes(y = reorder(.data[[order]], .data[[xval]]),
                 x = .data[[xval]],
                 fill = .data[[xval]])) +
      geom_col(show.legend = FALSE,
               color = "grey50") +
      scale_fill_gradient(
         low  = "white",
         high = "#2166AC"   # deep blue
      ) +
      geom_vline(
         xintercept = x_threshold,
         linetype = "dashed",
         linewidth = linewidth,
         alpha = alpha_line
      ) +
      labs(y = ylabel,
           x = xlabel) +
      ggpubr::theme_classic2() +
      theme(axis.title = element_text(size = 22),
            axis.ticks = element_blank())
   
   return(pl)
}

library(Seurat)
run_seurat_workflow <- function(
      counts,
      metadata = NULL,
      hvg = NULL,
      nfeatures = 2000,
      npcs = 30,
      do.center = TRUE,
      do.scale = TRUE,
      sketch = NULL,
      seed = 22
) {
   
   set.seed(seed)
   seu <- CreateSeuratObject(
      counts = counts,
      meta.data = metadata
   )
   seu <- NormalizeData(
      seu,
      normalization.method = "LogNormalize",
      scale.factor = 10000,
      verbose = FALSE
   )
   
   if(!is.null(sketch)){
      seu <- SketchData(
         object = seu,
         ncells = sketch,
         method = 'LeverageScore',
         sketched.assay = 'sketch', 
         seed = seed
         )
      
      DefaultAssay(seu) <- 'sketch'
   }
   
   if(is.null(hvg)){
      seu <- FindVariableFeatures(
         seu,
         selection.method = "vst",
         nfeatures = nfeatures,
         verbose = FALSE
      )
      
   } else {
      VariableFeatures(seu) <- hvg
   }
   
   seu <- ScaleData(
      seu,
      features = VariableFeatures(seu),
      do.center = do.center,
      do.scale = do.scale,
      verbose = FALSE
   )
   seu <- RunPCA(
      seu,
      features = VariableFeatures(seu),
      npcs = npcs,
      seed.use = 22,
      verbose = FALSE
   )
   seu <- RunUMAP(
      seu,
      dims = 1:npcs,
      seed.use = seed,
      verbose = FALSE
   )
   
   return(seu)
}


get_umap_embeddings <- function(scTypeEval,
                                assay = "single-cell",
                                seed = 22,
                                ...) {
   embeds <- scTypeEval@reductions[[assay]]@embeddings
   
   um <- uwot::umap(t(embeds), seed = seed, ...)
   colnames(um) <- c("umap_1", "umap_2")
   
   md <- scTypeEval@metadata
   md <- md[rownames(um),]
   
   if(!identical(rownames(md), rownames(um))){
      stop("Not identical rownames between umap embeddings and metadata")
   }
   df <- cbind(um, md)

   return(df)
}


plot_umap_embeddings <- function(df,
                                 color.by = "celltype",
                                 label = TRUE,
                                 label.by = color.by,
                                 colors = NULL,
                                 dot.size = 0.5,
                                 label.size = 8,
                                 show.legend = FALSE,
                                 keep.square = TRUE) {

   pl <- ggplot2::ggplot(df) +
      ggplot2::geom_point(
         ggplot2::aes(x = umap_1, y = umap_2, color = .data[[color.by]]),
         alpha = 0.8,
         size = dot.size,
         show.legend = show.legend
      ) +
      ggpubr::theme_classic2() +
      ggplot2::theme(
         legend.text = ggplot2::element_text(size = 12, face = "bold"),
         axis.title = ggplot2::element_text(size = 22),
         axis.text = ggplot2::element_blank(),
         axis.ticks = ggplot2::element_blank(),
         axis.line = ggplot2::element_blank()
      )
   
   if (keep.square) {
      pl <- pl + ggplot2::coord_fixed() + ggplot2::theme(aspect.ratio = 1)
   }
   
   if (!is.null(colors)) {
      pl <- pl + ggplot2::scale_color_manual(values = colors)
   }
   
   if (label) {
      cnt <- df |>
         dplyr::group_by(.data[[label.by]]) |>
         dplyr::summarize(umap_1 = mean(umap_1), umap_2 = mean(umap_2))
      pl <- pl +
         ggrepel::geom_label_repel(
            data = cnt,
            ggplot2::aes(x = umap_1, y = umap_2,
                         label = .data[[label.by]],
                         color = .data[[label.by]]),
            alpha = 0.9,
            size = label.size,
            max.overlaps = Inf,
            show.legend = FALSE
         )
   }
   return(pl)
}


get_pseudobulk <- function(scTypeEval,
                           filter = NULL,
                           genes = NULL,
                           order = NULL,
                           ident = "celltype",
                           sample = "sample",
                           vars = NULL){
   
   if(is.null(filter)){
      md <- scTypeEval@metadata
      
   } else {
      md <- scTypeEval@metadata %>% 
         filter(.data[[names(filter)]] == filter)
   }
   
   mat <- scTypeEval@data$pseudobulk@matrix
   cts <- scTypeEval@data$pseudobulk@ident[[ident]] %>%
      unique()
   
   if(!is.null(genes)){
      g <- intersect(genes, rownames(mat))
      mat <- mat[g,]
   }
   
   
   # get metadata
   ps_md <- md %>% 
      mutate(!!ident := scTypeEval:::purge_label(.data[[ident]]),
             !!sample := scTypeEval:::purge_label(.data[[sample]]),
            sample_id  = paste(.data[[sample]], .data[[ident]], sep = "_")
             ) %>% 
      filter(.data[[ident]] %in% cts,
             sample_id %in% colnames(mat))
   selected_vars <- c("sample_id", ident, sample, vars)
   ps_md <- ps_md[,selected_vars] %>% 
      distinct(sample_id, .keep_all = T)
   rownames(ps_md) <- ps_md$sample_id
   ps_md <- ps_md %>% 
      select(-sample_id)
   
   if(!is.null(order)){
      ps_md[[ident]] <- factor(ps_md[[ident]],
                               levels = order)
      ps_md <- ps_md %>% 
         arrange(.data[[ident]],
                 .data[[sample]])
   }
   
   mat <- mat[,rownames(ps_md)]
   ps_md <- ps_md %>% 
      select(-.data[[sample]])
   
   ret <- list(matrix = mat,
               metadata = ps_md)
}

library(pheatmap)
do_heatmap <- function(ps,
                       do_scale = TRUE,
                       cap = 4,
                       cluster_rows = TRUE,
                       cluster_cols = FALSE){
   
   if(do_scale){
      mat_scaled <- t(scale(t(ps$matrix)))
      mat_scaled <- pmax(pmin(mat_scaled, cap), -cap)
   } else {
      mat_scaled <- ps$matrix
   }
   
   gaps <- table(ps$metadata[[ident]]) %>% as.vector()
   gaps <- cumsum(gaps)
   
   ph <- pheatmap(
      mat = mat_scaled,
      cluster_cols = cluster_cols,
      cluster_rows = cluster_rows,
      annotation_col = ps$metadata,
      show_colnames = FALSE,
      gaps_col = gaps
   )
   
   return(ph)
   
}


sample_boxplot_expression <- function(seu,
                                      signature_pattern = "_score",
                                      ident = "celltype",
                                      ident_sel = "CD4.Tstr",
                                      sample = "sample") {
   
   df <- seu@meta.data %>% 
      filter(.data[[ident]] == ident_sel)
   keep <- grep(signature_pattern, names(df), value = T)
   keep <- c(sample, keep)
   
   df <- df[,keep]
   
   df <- df %>% 
      pivot_longer(-c(.data[[sample]]),
                   names_to = "signature",
                   values_to = "score") %>% 
      mutate(signature = factor(
         gsub("_score", "", signature),
         levels = sort(unique(gsub("_score", "", signature)))
      ))
   
   pl <- df %>% 
      ggplot(aes(.data[[sample]], score,
                 fill = signature)) +
      geom_boxplot(show.legend = F) +
      ylim(c(0,1)) +
      facet_wrap(~ signature,
                 ncol = 1) +
      labs(y = "Cell type gene signature score",
           x = sample,
           title = ident_sel) +
      ggpubr::theme_classic2() +
      theme(axis.text.x = element_text(angle = 45,
                                       hjust = 1,
                                       vjust = 1))
   
   
}

library(ProjecTILs)
self_project_pt <- function(
      seu,
      heldout = "CD4.Tstr",
      ident = "celltype",
      sample_split = 0.5,
      sample_id = "patient",
      npcs = 30,
      nfeatures = 2000,
      seed = 22
) {
   
   stopifnot(sample_split > 0, sample_split <= 1)
   
   meta <- seu[[]]
   
   # Validate metadata
   if (!all(c(ident, sample_id) %in% colnames(meta))) {
      stop("ident or sample_id is not present in Seurat metadata.")
   }
   
   if (!heldout %in% meta[[ident]]) {
      stop("Held-out cell type not found in metadata.")
   }
   
   # Select patients containing at least one held-out cell
   target_samples <- unique(
      as.character(meta[[sample_id]][
         meta[[ident]] == heldout &
            !is.na(meta[[sample_id]])
      ])
   )
   
   if (length(target_samples) < 2) {
      stop("At least two patients containing the held-out cell type are needed.")
   }
   
   # Hold out a subset of target-positive patients
   set.seed(seed)
   n_query <- min(
      length(target_samples) - 1L,
      ceiling(length(target_samples) * sample_split)
   )
   
   query_samples <- sample(target_samples, n_query)
   
   query_cells <- rownames(meta)[
      meta[[ident]] == heldout &
         as.character(meta[[sample_id]]) %in% query_samples
   ]
   
   # Exclude all cells from query patients from the reference
   ref_cells <- rownames(meta)[
      meta[[ident]] != heldout &
         !as.character(meta[[sample_id]]) %in% query_samples
   ]
   
   if (length(query_cells) == 0 || length(ref_cells) == 0) {
      stop("The query or reference contains no cells.")
   }
   
   # Create query from raw counts
   query <- subset(seu, cells = query_cells)
   
   query <- CreateSeuratObject(
      counts = GetAssayData(query, assay = "RNA", layer = "counts"),
      meta.data = query[[]]
   )
   
   # Build reference from remaining patients and cell types
   ref0 <- subset(seu, cells = ref_cells)
   
   ref <- run_seurat_workflow(
      counts = GetAssayData(ref0, assay = "RNA", layer = "counts"),
      metadata = ref0[[]],
      npcs = npcs
   )
   
   # Create ProjecTILs reference
   ref_manual <- ProjecTILs::make.reference(
      ref = ref,
      assay = "RNA",
      ndim = npcs,
      seed = seed,
      recalculate.umap = TRUE,
      nfeatures = nfeatures,
      annotation.column = ident
   )
   
   # Project held-out cells onto the reference
   proj <- ProjecTILs::Run.ProjecTILs(
      query = query,
      ref = ref_manual,
      filter.cells = FALSE,
      ndim = npcs,
      ncores = 1,
      progressbar = FALSE,
      fast.umap.predict = FALSE
   )
   
   # UMAP visualization
   p_projection <- ProjecTILs::plot.projection(
      ref = ref_manual,
      query = proj,
      linesize = 0.2,
      pointsize = 1.2
   ) +
      ggplot2::ggtitle(paste("Held out:", heldout))
   
   # Which reference annotations receive the query cells?
   p_composition <- ProjecTILs::plot.statepred.composition(
      ref = ref_manual,
      query = proj
   ) +
      ggplot2::ggtitle(paste("Predicted states:", heldout))
   
   return(list(
      projection = p_projection,
      composition = p_composition,
      projected_query = proj,
      reference = ref_manual,
      query_samples = query_samples
   ))
}


library(SingleR)
library(SingleCellExperiment)

reannotate_singleR <- function(seu,
                               heldout = "CD4.Tstr",
                               ident = "celltype",
                               ncores = 1) {
   
   # Query: held-out cell type in held-out samples
   query_cells <- colnames(seu)[
      seu@meta.data[[ident]] == heldout]
   
   # Reference: all other cell types in remaining samples
   ref_cells <- colnames(seu)[
      seu@meta.data[[ident]] != heldout]
   
   if (length(query_cells) == 0 || length(ref_cells) == 0) {
      stop("No query or reference cells available. Check the split and metadata.")
   }
   
   # Extract raw counts
   counts <- GetAssayData(seu, assay = "RNA", layer = "counts")
   
   # Create SingleCellExperiment objects
   query <- SingleCellExperiment(
      assays = list(counts = counts[, query_cells, drop = FALSE])
   )
   query <- scuttle::logNormCounts(query)
   
   ref <- SingleCellExperiment(
      assays = list(counts = counts[, ref_cells, drop = FALSE])
   )
   ref <- scuttle::logNormCounts(ref)
   
   # Reference labels
   ref_labels <- seu@meta.data[ref_cells, ident, drop = TRUE]
   ref_labels <- as.character(ref_labels)
   
   bparam <- MulticoreParam(workers = ncores,
                            progressbar = F)
   
   # Run SingleR
   pred <- SingleR(
      test = query,
      ref = ref,
      labels = ref_labels,
      de.method = "classic",
      BPPARAM = bparam
   )
   
   # Return predictions and evaluation data
   return(as.data.frame(pred))
}



dotplot_seurat_style <- function(
      object,
      features,
      ident_col = "annotation",
      split_by = "tissue",
      subset_col = NULL,
      subset_values = NULL,
      assay = "RNA",
      layer = "data",
      scale_expression = TRUE,
      low_color = "gray80",
      high_color = "red3",
      dot_scale = 7
) {
   stopifnot(inherits(object, "Seurat"))
   
   md <- object@meta.data
   if (!ident_col %in% colnames(md)) {
      stop(paste0("ident_col not found in metadata: ", ident_col))
   }
   
   if (!is.null(split_by) && !split_by %in% colnames(md)) {
      stop(paste0("split_by not found in metadata: ", split_by))
   }
   
   if (!is.null(subset_col) && !subset_col %in% colnames(md)) {
      stop(paste0("subset_col not found in metadata: ", subset_col))
   }
   
   if (is.null(subset_col) != is.null(subset_values)) {
      stop("subset_col and subset_values must be provided together.")
   }
   
   if (!is.null(subset_col)) {
      keep_cells <- rownames(md)[md[[subset_col]] %in% subset_values]
      if (length(keep_cells) == 0) {
         stop("No cells found after applying subset_col/subset_values.")
      }
      object <- subset(object, cells = keep_cells)
      md <- object@meta.data
   }
   
   if (!is.null(assay)) {
      DefaultAssay(object) <- assay
   }
   
   feat_found <- features[features %in% rownames(object)]
   if (length(feat_found) == 0) {
      stop("None of the requested features were found in the object.")
   }
   
   if (length(feat_found) < length(features)) {
      missing_feats <- setdiff(features, feat_found)
      message("Skipping missing features: ", paste(missing_feats, collapse = ", "))
   }
   
   vars_to_fetch <- c(ident_col, split_by, feat_found) %>%
      unique() %>%
      .[!is.na(.)]
   
   expr_df <- Seurat::FetchData(object = object, vars = vars_to_fetch, layer = layer)
   
   plot_df <- expr_df %>%
      tibble::rownames_to_column("cell") %>%
      tidyr::pivot_longer(
         cols = tidyr::all_of(feat_found),
         names_to = "feature",
         values_to = "expr"
      ) %>%
      dplyr::mutate(
         ident = .data[[ident_col]],
         split = if (!is.null(split_by)) .data[[split_by]] else "all"
      ) %>%
      dplyr::group_by(ident, split, feature) %>%
      dplyr::summarise(
         pct_exp = mean(expr > 0) * 100,
         avg_exp = mean(expr),
         .groups = "drop"
      ) %>%
      dplyr::group_by(feature) %>%
      dplyr::mutate(
         avg_exp_scaled = as.numeric(scale(avg_exp)),
         avg_exp_scaled = ifelse(is.na(avg_exp_scaled), 0, avg_exp_scaled)
      ) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(
         feature = factor(feature, levels = rev(feat_found)),
         ident = droplevels(factor(ident)),
         split = droplevels(factor(split)),
         color_value = if (scale_expression) avg_exp_scaled else avg_exp
      )
   
   # Rebuild ident levels from observed values only to avoid spacing artifacts.
   ident_levels <- unique(as.character(plot_df$ident))
   plot_df <- plot_df %>%
      dplyr::mutate(
         ident = factor(as.character(ident), levels = ident_levels),
         ident_idx = as.numeric(ident)
      )
   
   if (!is.null(split_by)) {
      n_split <- nlevels(plot_df$split)
      split_levels <- levels(plot_df$split)
      split_offsets <- seq_len(n_split) - (n_split + 1) / 2
      names(split_offsets) <- split_levels
      
      plot_df <- plot_df %>%
         dplyr::mutate(
            split_idx = as.numeric(split),
            x_pos = ident_idx + split_offsets[as.character(split)] * 0.28
         )
      
      split_breaks <- as.vector(
         sapply(seq_len(nlevels(plot_df$ident)), function(i) i + split_offsets * 0.28)
      )
      split_labels <- rep(split_levels, times = nlevels(plot_df$ident))
   } else {
      plot_df <- plot_df %>%
         dplyr::mutate(x_pos = ident_idx)
   }
   
   p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = x_pos, y = feature)) +
      ggplot2::geom_point(
         ggplot2::aes(size = pct_exp, color = color_value)
      ) +
      ggplot2::scale_size(range = c(0, dot_scale), name = "% Percent\nExpressed") +
      ggplot2::scale_color_gradient(
         low = low_color,
         high = high_color,
         name = if (scale_expression) "Average\nExpression\n(scaled)" else "Average\nExpression"
      ) +
      ggplot2::labs(x = NULL, y = NULL) +
      ggplot2::theme_bw() +
      ggplot2::theme(
         axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 18),
         panel.grid.minor = ggplot2::element_blank()
      )
   
   if (!is.null(split_by)) {
      p <- p +
         ggplot2::geom_vline(
            xintercept = seq_len(nlevels(plot_df$ident)) + 0.5,
            color = "gray90",
            linewidth = 0.3
         ) +
         ggplot2::scale_x_continuous(
            breaks = split_breaks,
            labels = split_labels,
            sec.axis = ggplot2::dup_axis(
               breaks = seq_len(nlevels(plot_df$ident)),
               labels = levels(plot_df$ident),
               name = ident_col
            )
         ) +
         ggplot2::labs(x = split_by) +
         ggplot2::theme(
            axis.title.x.top = ggplot2::element_text(margin = ggplot2::margin(b = 6)),
            axis.text.x.top = ggplot2::element_text(angle = 0, hjust = 0.5),
            axis.text.x = ggplot2::element_text(angle = 45,
                                                hjust = 1,
                                                size = 18)
         )
   } else {
      p <- p +
         ggplot2::scale_x_continuous(
            limits = c(0.5, nlevels(plot_df$ident) + 0.5),
            breaks = seq_len(nlevels(plot_df$ident)),
            labels = levels(plot_df$ident),
            expand = ggplot2::expansion(mult = c(0, 0))
         ) +
         ggplot2::labs(x = ident_col)
   }
   
   return(p)
}


remove_batch_effect <- function(
      mat,
      batch,
      covariates = NULL
) {
   stopifnot(ncol(mat) == length(batch))
   
   design <- NULL
   
   if (!is.null(covariates)) {
      design <- model.matrix(~ ., data = covariates)
   }
   
   limma::removeBatchEffect(
      x = as.matrix(mat),
      batch = factor(batch),
      design = design
   )
}


process_metadata_pseudobulk <- function(md,
                                        ident = "cell_type",
                                        sample = "sample",
                                        covariates = NULL){
   df <- md %>% 
      mutate(!!sample := scTypeEval:::purge_label(.data[[sample]]),
             !!ident := scTypeEval:::purge_label(.data[[ident]]),
             id = paste(.data[[sample]], .data[[ident]], sep = "_")
      ) %>% 
      .[,c("id", covariates)] %>% 
      distinct() %>% 
      tibble::remove_rownames() %>% 
      tibble::column_to_rownames("id")
   
   return(df)
}
