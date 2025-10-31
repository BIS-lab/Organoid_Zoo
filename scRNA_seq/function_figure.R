library(dplyr)
library(tibble)
library(pheatmap)
library(RColorBrewer)
library(ggthemes) 
library(paletteer)

plot_featureplot_v2 <- function(obj, gene, max){
  DefaultAssay(obj) <- "RNA"
  ht_custom_col <- rev(paletteer_c("grDevices::RdPu", 30))
  p11 <- FeaturePlot(obj, features = gene, pt.size = 1, min.cutoff = 0.0, max.cutoff = max, order = TRUE)
  ggplot(data = p11$data, aes_string(x = 'UMAP_1', y = 'UMAP_2', fill = gene)) +    
    geom_point(shape = 21, stroke = 0.3, color = "black", alpha = 1, size = 2) +
    scale_fill_gradientn(colours = ht_custom_col, limits=c(0.0, max)) +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
    labs(fill = "Stem cell \nsignature score")
}

plot_featureplot_v3 <- function(obj, gene, min, max){
  DefaultAssay(obj) <- "RNA"
  ht_custom_col <- rev(paletteer_c("grDevices::Spectral", 30))
  p11 <- FeaturePlot(obj, features = gene, pt.size = 1, min.cutoff = min, max.cutoff = max, order = TRUE)
  ggplot(data = p11$data, aes_string(x = 'UMAP_1', y = 'UMAP_2', fill = gene)) +    
    geom_point(shape = 21, stroke = 0.3, color = "black", alpha = 1, size = 2) +
    scale_fill_gradientn(colours = ht_custom_col, limits=c(min, max)) +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
}

plot_featureplot_v4 <- function(obj, gene, max){
  DefaultAssay(obj) <- "RNA"
  base_cols <- rev(paletteer_c("grDevices::RdPu", 30))
  low_col   <- base_cols[1]          
  high_col  <- "#B82647"            

  p11 <- FeaturePlot(
    obj, features = gene, pt.size = 1,
    min.cutoff = 0.0, max.cutoff = max, order = TRUE
  )
  ggplot(p11$data, aes_string(x = "UMAP_1", y = "UMAP_2", fill = gene)) +
    geom_point(shape = 21, stroke = 0.01, color = "black", alpha = 1, size = 2) +
    scale_fill_gradient(limits = c(0.0, max), low = low_col, high = high_col) +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5)
    ) +
    labs(fill = "CBC \nsignature score")
}

plot_featureplot_v5 <- function(obj, gene, max){
  DefaultAssay(obj) <- "RNA"
  base_cols <- rev(paletteer_c("grDevices::PuBu", 30))
  low_col   <- base_cols[1]         
  high_col  <- "#0B6DB7"             
  p11 <- FeaturePlot(
    obj, features = gene, pt.size = 1,
    min.cutoff = 0.0, max.cutoff = max, order = TRUE
  )
  ggplot(p11$data, aes_string(x = "UMAP_1", y = "UMAP_2", fill = gene)) +
    geom_point(shape = 21, stroke = 0.01, color = "black", alpha = 1, size = 2) +
    scale_fill_gradient(limits = c(0.0, max), low = low_col, high = high_col) +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5)
    ) +
    labs(fill = "revSC \nsignature score")
}