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
      nfeatures = 2000,
      npcs = 30,
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
   seu <- FindVariableFeatures(
      seu,
      selection.method = "vst",
      nfeatures = nfeatures,
      verbose = FALSE
   )
   seu <- ScaleData(
      seu,
      features = VariableFeatures(seu),
      verbose = FALSE
   )
   seu <- RunPCA(
      seu,
      features = VariableFeatures(seu),
      npcs = npcs,
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
