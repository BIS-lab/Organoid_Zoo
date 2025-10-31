library(tibble)
library(dplyr)
library(ggplot2)
library(RColorBrewer)

averageMat <- function(seurat, Genes, group.by = "seurat_clusters", scale=TRUE){
  
  selected_geneids <- Genes
  
  # Select only valid genes
  valid_genes <- selected_geneids[selected_geneids %in% rownames(seurat@assays$RNA@data)]
  
  if(length(valid_genes) == 0){
    stop("No valid genes found in the matrix.")
  }
  
  data <- as.matrix(seurat@assays$RNA@data[valid_genes, ])
  if(sum(scale) == 1){scale.data <- t(scale(t(data)))}else{scale.data <- data}
  
  df <- data.frame(row.names = colnames(seurat), group=seurat@meta.data[,group.by], t(scale.data))
  df <- df %>%
    group_by(group) %>%
    summarise_at(vars(all_of(valid_genes)), list(mean))
  df <- tibble::column_to_rownames(df, var="group")
  colnames(df) <- valid_genes
  mat <- as.matrix(df)
  
  if(sum(scale) == 1){
    mat[mat > 1.5] <- 1.5
    mat[mat < -1.5] <- -1.5
  }
  
  return(mat)
}

propMat <- function(seurat, Genes, group.by = "seurat_clusters"){
  
  selected_geneids <- Genes
  
  # Select only valid genes
  valid_genes <- selected_geneids[selected_geneids %in% rownames(seurat@assays$RNA@data)]
  
  if(length(valid_genes) == 0){
    stop("No valid genes found in the matrix.")
  }
  
  data <- as.matrix(seurat@assays$RNA@data[valid_genes, ])
  propdata <- data > 0
  propdata <- propdata * 1
  
  df <- data.frame(row.names = colnames(seurat), group=seurat@meta.data[,group.by], t(propdata))
  df <- df %>%
    group_by(group) %>%
    summarise_at(vars(all_of(valid_genes)), list(mean))
  df <- tibble::column_to_rownames(df, var="group")
  colnames(df) <- valid_genes
  prop.mat <- as.matrix(df)
  
  return(prop.mat)
}

average_dotplot <- function(seurat, Genes, group.by = "seurat_clusters", cl.order=NULL, cells=NULL, range=c(0,6)){
  
  if(sum(is.null(cells)) == 0){
    seurat <- subset(seurat, cells=cells)
  }
  
  if(sum(is.null(cl.order)) == 0){
    cl.level <- cl.order
    seurat@meta.data[, group.by] <- factor(seurat@meta.data[, group.by], levels = cl.order)
  }else{
    if(sum(is.null(levels(seurat@meta.data[, group.by]))) == 0){
      cl.level <- levels(seurat@meta.data[, group.by])
    }else{
      cl.level <- unique(seurat@meta.data[, group.by])
    }
  }
  
  # Remove duplicate genes
  unique_genes <- unique(Genes)
  
  # Select only valid genes
  valid_genes <- unique_genes[unique_genes %in% rownames(seurat@assays$RNA@data)]
  
  # Non-existent genes are treated as NA
  missing_genes <- unique_genes[!unique_genes %in% rownames(seurat@assays$RNA@data)]
  
  if(length(missing_genes) > 0){
    warning("The following genes are missing from the data: ", paste(missing_genes, collapse = ", "))
  }
  
  # Maintain the original gene order
  all_genes <- unique_genes
  
  mat <- averageMat(seurat, valid_genes, group.by)
  prop.mat <- propMat(seurat, valid_genes, group.by)
  
  # Check existing column names in `mat` and `prop.mat`
  existing_genes <- colnames(mat)
  
  # Expand data to include all genes (original order)
  mat <- mat[, existing_genes, drop = FALSE] # Existing valid_genes
  prop.mat <- prop.mat[, existing_genes, drop = FALSE]
  
  # Add missing genes as NA
  for (gene in missing_genes) {
    mat <- cbind(mat, setNames(data.frame(rep(NA, nrow(mat))), gene))
    prop.mat <- cbind(prop.mat, setNames(data.frame(rep(NA, nrow(prop.mat))), gene))
  }
  
  # Sort the column order in original order (all_genes)
  mat <- mat[, all_genes, drop = FALSE]
  prop.mat <- prop.mat[, all_genes, drop = FALSE]
  
  dot.df <- data.frame(row.names = cl.level)
  dot.df <- cbind.data.frame(dot.df, mat)
  
  data.df <- data.frame(genes = rep(all_genes, length(cl.level)))
  data.df$celltype <- rep(cl.level, each=length(all_genes))
  data.df$expr <- as.vector(t(mat))
  data.df$prop <- as.vector(t(prop.mat))
  
  min.value <- min(data.df$expr, na.rm = TRUE)
  max.value <- max(data.df$expr, na.rm = TRUE)
  abs.value <- max(abs(min.value), max.value, na.rm = TRUE)
  
  data.df$celltype <- as.character(data.df$celltype)
  
  gd <- ggplot(data.df) + 
    geom_point(aes(x=genes, y=celltype, color=expr, size=prop)) +
    scale_x_discrete(limits = all_genes) + 
    scale_y_discrete(limits = as.character(cl.level)) +
    scale_color_gradientn(colors = rev(c(brewer.pal(11, "RdBu")[1:5], rep(brewer.pal(11, "RdBu")[6], 2), brewer.pal(11, "RdBu")[7:11])), limits = c(-abs.value, abs.value), na.value = "white") + 
    scale_radius(range = range) +
    theme_classic() +
    theme(text = element_text(size=20)) +
    ylab("") + xlab("") + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  
  return(gd)
}
