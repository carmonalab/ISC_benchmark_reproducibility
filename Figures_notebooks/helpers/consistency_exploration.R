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
