library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(SeuratDisk)
library(tibble)
library(pheatmap)
library(RColorBrewer)
library(ggthemes) 
library(paletteer)
library(matrixStats)
library(purrr)

source("function_figure.R")

### Load singlet objects for all species
mms <- LoadH5Seurat("../data/mouse_Ms_singlet.h5Seurat")
mh <- LoadH5Seurat("../data/mouse_H_singlet.h5Seurat")
bovine <- LoadH5Seurat("../data/bovine_singlet.h5Seurat")
ferret <- LoadH5Seurat("../data/ferret_singlet.h5Seurat")
spiny <- LoadH5Seurat("../data/spiny_singlet.h5Seurat")
rat <- LoadH5Seurat("../data/rat_singlet.h5Seurat")
chick <- LoadH5Seurat("../data/chick_singlet.h5Seurat")
beagle <- LoadH5Seurat("../data/beagle_singlet.h5Seurat")
hamster <- LoadH5Seurat("../data/hamster_singlet.h5Seurat")
hmouse <- LoadH5Seurat("../data/hmouse_singlet.h5Seurat")
human <- LoadH5Seurat("../data/human_singlet.h5Seurat")
marmoset <- LoadH5Seurat("../data/marmoset_singlet.h5Seurat")
sakhalin <- LoadH5Seurat("../data/sakhalin_singlet.h5Seurat")
steppe <- LoadH5Seurat("../data/steppe_singlet.h5Seurat")
pig <- LoadH5Seurat("../data/pig_singlet.h5Seurat")
goat <- LoadH5Seurat("../data/goat_singlet.h5Seurat")
hypsugo <- LoadH5Seurat("../data/hypsugo_bat.h5Seurat")


# -----------------------------------------------------------------------------
# Utility and plotting functions for signature scoring & cross-species analysis
# -----------------------------------------------------------------------------


# --- Species registry (Seurat objects + ortholog mapping info) ---------------
#  Load singlet objects for all species
objs <- list(
  mms       = mms,
  mh        = mh,
  bovine    = bovine,
  ferret    = ferret,
  spiny     = spiny,
  rat       = rat,
  chick     = chick,
  beagle    = beagle,
  hamster   = hamster,
  hmouse    = hmouse,
  human     = human,
  marmoset  = marmoset,
  sakhalin  = sakhalin,
  steppe    = steppe,
  pig       = pig,
  goat      = goat,
  hypsugo   = hypsugo
)

# Load signature gene sets from Extended Data Table 2 
# CBC stem-cell, revival stem-cell, and secretory signatures.
signature.raw <- read.csv("../data/extended_data_table2.csv", header = TRUE)
stem.sig <- na.omit(signature.raw$CBC_stem_cell_signature_geneset)
rev.sig  <- na.omit(signature.raw$revival_stem_cell_signature_geneset)
sec.sig  <- na.omit(signature.raw$secretory_signature_geneset)

# Load ortholog mapping files
ferret.df <- read.csv(file = "../data/mms_mpfuro_ortholog.csv")
ferret.df.one2one <- ferret.df[ferret.df$mpfuro_homolog_orthology_type == "ortholog_one2one",]

bovine.df <- read.csv(file = "/vol/data/mms_btaurus_ortholog.csv")
bovine.df.one2one <- bovine.df[bovine.df$btaurus_homolog_orthology_type == "ortholog_one2one",]

rat.df <- read.csv(file = "/vol/data/mms_rnorvegicus_ortholog.csv")
rat.df.one2one <- rat.df[rat.df$rnorvegicus_homolog_orthology_type == "ortholog_one2one",]

chick.df <- read.csv(file = "/vol/data/mms_ggallus_ortholog.csv")
chick.df.one2one <- chick.df[chick.df$ggallus_homolog_orthology_type == "ortholog_one2one",]

beagle.df <- read.csv(file = "/vol/data/mms_clfamiliaris_ortholog.csv")
beagle.df.one2one <- beagle.df[beagle.df$clfamiliaris_homolog_orthology_type == "ortholog_one2one",]

hamster.df <- read.csv(file = "/vol/data/mms_mauratus_ortholog.csv")
hamster.df.one2one <- hamster.df[hamster.df$mauratus_homolog_orthology_type == "ortholog_one2one",]

human.df <- read.csv(file = "/vol/data/mms_hsapiens_ortholog.csv")
human.df.one2one <- human.df[human.df$hsapiens_homolog_orthology_type == "ortholog_one2one",]

marmoset.df <- read.csv(file = "/vol/data/mms_cjacchus_ortholog.csv")
marmoset.df.one2one <- marmoset.df[marmoset.df$cjacchus_homolog_orthology_type == "ortholog_one2one",]

pig.df <- read.csv(file = "/vol/data/mms_sscrofa_ortholog.csv")
pig.df.one2one <- pig.df[pig.df$sscrofa_homolog_orthology_type == "ortholog_one2one",]

goat.df <- read.csv(file = "/vol/data/mms_chircus_ortholog.csv")
goat.df.one2one <- goat.df[goat.df$chircus_homolog_orthology_type == "ortholog_one2one",]

# Load ortholog mapping files; see Materials and Methods for details
spiny.df    <- read.table(file = "/vol/spiny/spiny_mouse_ortholog_final.txt"); rownames(spiny.df) <- spiny.df$V1
hmouse.df   <- read.csv(file = "/vol/data/hmouse_mouse_ortholog_final.csv"); rownames(hmouse.df) <- hmouse.df$V1
sakhalin.df <- read.csv(file = "/vol/data/sakhalin_mouse_ortholog_final.csv"); rownames(sakhalin.df) <- sakhalin.df$V1
steppe.df   <- read.csv(file = "/vol/data/steppe_mouse_ortholog_final.csv");   rownames(steppe.df) <- steppe.df$V1
hypsugo.df  <- read.csv(file = "/vol/data/hypsugo_mouse_ortholog_final.csv");  rownames(hypsugo.df) <- hypsugo.df$V1

# Species → unified ortholog registry (from mouse to target species)
orth <- list(
  mms      = NULL,   # mouse: no mapping required
  mh       = NULL,   # mouse H: no mapping required
  rat      = list(df = rat.df.one2one,      from_col = "external_gene_name", to_col = "rnorvegicus_homolog_associated_gene_name"),
  spiny    = list(df = spiny.df,            from_col = "V1",                 to_col = "V2"),
  ferret   = list(df = ferret.df.one2one,   from_col = "external_gene_name", to_col = "mpfuro_homolog_associated_gene_name"),
  bovine   = list(df = bovine.df.one2one,   from_col = "external_gene_name", to_col = "btaurus_homolog_associated_gene_name"),
  chick    = list(df = chick.df.one2one,    from_col = "external_gene_name", to_col = "ggallus_homolog_associated_gene_name"),
  beagle   = list(df = beagle.df.one2one,   from_col = "external_gene_name", to_col = "clfamiliaris_homolog_associated_gene_name"),
  hamster  = list(df = hamster.df.one2one,  from_col = "external_gene_name", to_col = "mauratus_homolog_associated_gene_name"),
  human    = list(df = human.df.one2one,    from_col = "external_gene_name", to_col = "hsapiens_homolog_associated_gene_name"),
  hmouse   = list(df = hmouse.df,           from_col = "V1",                 to_col = "V2"),
  marmoset = list(df = marmoset.df.one2one, from_col = "external_gene_name", to_col = "cjacchus_homolog_associated_gene_name"),
  pig      = list(df = pig.df.one2one,      from_col = "external_gene_name", to_col = "sscrofa_homolog_associated_gene_name"),
  goat     = list(df = goat.df.one2one,     from_col = "external_gene_name", to_col = "chircus_homolog_associated_gene_name"),
  sakhalin = list(df = sakhalin.df,         from_col = "V1",                 to_col = "V2"),
  steppe   = list(df = steppe.df,           from_col = "V1",                 to_col = "V2"),
  hypsugo  = list(df = hypsugo.df,          from_col = "V1",                 to_col = "V2")
)

# --- Core utility functions ---------------------------------------------------

# Map mouse gene set to a target species based on ortholog table
map_genes_to_species <- function(mouse_genes, species, orth_info, obj) {
  if (is.null(orth_info)) {
    return(intersect(mouse_genes, rownames(obj)))
  }
  df <- orth_info$df
  from_col <- orth_info$from_col
  to_col   <- orth_info$to_col
  if (!has_rows(df)) return(character(0))
  df2 <- df[, c(from_col, to_col)] %>% distinct() %>% na.omit()
  mapped <- df2 %>%
    filter(.data[[from_col]] %in% mouse_genes) %>%
    pull(.data[[to_col]]) %>% unique()
  intersect(mapped, rownames(obj))
}

# AddModuleScore wrapper
add_signature_score <- function(obj, genes, name = "stem", assay = "RNA") {
  if (length(genes) == 0) return(obj)
  AddModuleScore(obj, features = list(genes), name = name, assay = assay)
}

# Compute correlation between signature score and expression (Spearman by default)
cor_vs_score <- function(obj, score_col = "stem1", genes, method = "spearman", slot = "data") {
  if (!score_col %in% colnames(obj@meta.data)) stop("score_col not found in meta.data")
  if (length(genes) == 0) return(setNames(numeric(0), character(0)))
  mat <- GetAssayData(obj, slot = slot)
  genes <- intersect(genes, rownames(mat))
  if (length(genes) == 0) return(setNames(numeric(0), character(0)))

  ord <- order(-obj@meta.data[[score_col]])
  score_raw <- obj@meta.data[[score_col]][ord]
  # Scale score to [0,1] 
  score <- (score_raw - min(score_raw)) / (max(score_raw) - min(score_raw) + 1e-12)

  expr <- as.matrix(mat[genes, , drop = FALSE])[, ord, drop = FALSE]
  cors <- apply(expr, 1, function(x) suppressWarnings(cor(score, x, method = method)))
  as.numeric(cors) %>% setNames(rownames(expr))
}

# Bind species-wise correlation vectors into a single matrix in mouse-gene order
bind_cor_tables <- function(cor_list, mouse_ref_genes) {
  dfs <- imap(cor_list, function(v, nm) {
    tibble(gene = names(v), !!nm := as.numeric(v))
  }) %>%
    reduce(full_join, by = "gene") %>%
    filter(gene %in% mouse_ref_genes) %>%
    arrange(match(gene, mouse_ref_genes)) %>%
    column_to_rownames("gene")
  dfs
}

# Row ordering by (mean - sd) score
order_by_mean_minus_sd <- function(mat) {
  m <- rowMeans(mat, na.rm = TRUE)
  s <- rowSds(mat, na.rm = TRUE)
  ord <- order(m - s, decreasing = TRUE)
  mat[ord, , drop = FALSE]
}

order_rows_mean_minus_sd_with_NA <- function(x) {
  # coerce to data.frame, keep numerics only (row ops need numeric)
  x_df <- as.data.frame(x, check.names = FALSE)
  num_cols <- vapply(x_df, is.numeric, logical(1))
  if (!any(num_cols)) {
    # nothing numeric: return as-is
    return(x)
  }
  mat <- as.matrix(x_df[, num_cols, drop = FALSE])

  cc <- stats::complete.cases(mat)
  mat_noNA <- mat[cc, , drop = FALSE]
  mat_NA   <- mat[!cc, , drop = FALSE]

  # (mean - sd) for complete rows
  if (nrow(mat_noNA) > 0) {
    m <- rowMeans(mat_noNA)
    s <- matrixStats::rowSds(mat_noNA)
    score <- m - s
    ord_noNA <- order(score, decreasing = TRUE)
    r1 <- rownames(mat_noNA)[ord_noNA]
  } else {
    r1 <- character(0)
  }

  # NA rows: fewest NA first
  if (nrow(mat_NA) > 0) {
    na_counts <- rowSums(is.na(mat_NA))
    ord_NA <- order(na_counts, decreasing = FALSE)
    r2 <- rownames(mat_NA)[ord_NA]
  } else {
    r2 <- character(0)
  }

  final_order <- c(r1, r2)
  # Reindex original x in that order (preserve any non-numeric columns too)
  x[final_order, , drop = FALSE]
}

# Scatter of signature score and return cell order
signature_scatter <- function(obj, score_col = "stem1", ylim = c(-0.5, 1)) {
  meta <- obj@meta.data[order(obj@meta.data[[score_col]]), , drop = FALSE]
  df <- tibble(Cells = seq_len(nrow(meta)), Signature_score = meta[[score_col]])
  p <- ggplot(df, aes(Cells, Signature_score)) +
    geom_point(alpha = 0.7) +
    coord_cartesian(ylim = ylim) +
    theme(axis.text.x = element_blank())
  list(plot = p, cell_order = rownames(meta))
}

# Heatmap (z-scaled RNA@scale.data) using given gene & cell order
plot_heatmap_simple <- function(obj, gene_list, cell_order,
                                palette = colorRampPalette(c("#213156","#5790BF","white","#C96652","#5F1120"))(11),
                                breaks = seq(-3, 3, length.out = 12)) {
  sdata <- obj[["RNA"]]@scale.data
  missing <- setdiff(gene_list, rownames(sdata))
  sel <- sdata[rownames(sdata) %in% gene_list, , drop = FALSE]
  if (length(missing) > 0) {
    na_mat <- matrix(NA, nrow = length(missing), ncol = ncol(sdata),
                     dimnames = list(missing, colnames(sdata)))
    sel <- rbind(sel, na_mat)
  }
  sel <- sel[match(gene_list, rownames(sel)), cell_order, drop = FALSE]
  pheatmap(sel, color = palette, breaks = breaks,
           cluster_rows = FALSE, cluster_cols = FALSE,
           show_rownames = TRUE, show_colnames = FALSE)
}

prep_radar_data <- function(seurat_obj,
                            score_cols = c("stem1","secretory1","rev1"),
                            cluster_col = "celltype",
                            legend_order = names(my.cols)) {
  # 필수 컬럼 확인
  stopifnot(all(score_cols %in% colnames(seurat_obj@meta.data)))
  stopifnot(cluster_col %in% colnames(seurat_obj@meta.data))

  meta_pos <- seurat_obj@meta.data |>
    dplyr::select(all_of(c(cluster_col, score_cols))) |>
    mutate(across(all_of(score_cols), ~ ifelse(.x < 0, 0, .x)))

  df_tmp <- meta_pos |>
    group_by(.data[[cluster_col]]) |>
    summarise(across(all_of(score_cols), ~ mean(.x, na.rm = TRUE)), .groups = "drop") |>
    rename(group = !!cluster_col)

  if (nrow(df_tmp) == 0L) return(NULL)

  sig_stats <- df_tmp |>
    pivot_longer(all_of(score_cols), names_to = "Signature", values_to = "Score") |>
    group_by(Signature) |>
    summarise(min_val = min(Score, na.rm = TRUE),
              max_val = max(Score, na.rm = TRUE), .groups = "drop")

  df_scaled <- df_tmp |>
    pivot_longer(all_of(score_cols), names_to = "Signature", values_to = "Score") |>
    left_join(sig_stats, by = "Signature") |>
    mutate(Score_scaled = ifelse(max_val == min_val, 0, (Score - min_val)/(max_val - min_val))) |>
    select(group, Signature, Score_scaled) |>
    pivot_wider(names_from = Signature, values_from = Score_scaled)

  if (length(legend_order)) {
    lv <- c(intersect(legend_order, df_scaled$group),
            setdiff(df_scaled$group, legend_order))
    df_scaled$group <- factor(df_scaled$group, levels = lv)
    df_scaled <- df_scaled[order(df_scaled$group), , drop = FALSE]
  }

  as.data.frame(df_scaled, stringsAsFactors = FALSE)
}

plot_radar <- function(df_scaled, title, colors = my.cols) {
  if (is.null(df_scaled) || nrow(df_scaled) == 0L) return(NULL)

  grp <- as.character(df_scaled$group)


  cols_use <- colors[grp]
  if (any(is.na(cols_use))) {
    miss <- is.na(cols_use)
    auto_pal <- scales::hue_pal()(sum(miss))
    cols_use[miss] <- auto_pal
    names(cols_use) <- grp
  }

  ggradar(
    df_scaled,
    group.colours   = cols_use,
    grid.min        = 0,
    grid.mid        = 0.5,
    grid.max        = 1,
    values.radar    = c("0", "0.5", "1"),
    group.line.width = 1,
    group.point.size = 3,
    legend.position  = "right",
    axis.label.size  = 4
  ) +
    ggtitle(title) +
    scale_colour_manual(values = cols_use, drop = FALSE) +
    scale_fill_manual(values = alpha(cols_use, 0.2), drop = FALSE)
}



#################### stem cell signature gene score ###########################
mouse_stem <- unique(stem.sig)

# Map mouse stem signature to each species
mapped_stem <- imap(objs, ~ map_genes_to_species(mouse_stem, species = .y, orth_info = orth[[.y]], obj = .x))

# AddModuleScore for each species 
objs2 <- imap(objs, ~ add_signature_score(.x, mapped_stem[[.y]], name = "stem", assay = "RNA"))

# UMAP FeaturePlot with title including the number of signature genes
walk(names(objs2), function(sp) {
  n_genes <- length(mapped_stem[[sp]])
  p <- plot_featureplot_v4(objs2[[sp]], "stem1", 0.8) +
    labs(title = sprintf("Stem signature in %s (#: %d)", sp, n_genes)) +
    theme(plot.title = element_text(size = 20, face = "bold"))
  ggsave(p, filename = sprintf("/vol/output/2025final/stem_signature_%s_umap.pdf", sp),
         width = 5, height = 3, dpi = 300) #Extended Data Fig. 3d
})

# Calculate correlation (score vs. expression) per species
cors_stem <- imap(objs2, ~ cor_vs_score(.x, score_col = "stem1", genes = mapped_stem[[.y]], method = "spearman"))

# Combine species into two panels (as in original union/other species split)
cor_mat_A <- bind_cor_tables(cors_stem[c("mms","mh","rat","spiny","ferret","bovine","chick")], mouse_ref_genes = mouse_stem)
cor_mat_B <- bind_cor_tables(cors_stem[c("hmouse","hamster","human","marmoset","beagle","pig","goat","hypsugo","steppe","sakhalin")], mouse_ref_genes = mouse_stem)

# Heatmap parameters 
my_colors <- colorRampPalette(c("#213156","#5790BF","white","#C96652","#5F1120"))(101)
my_breaks <- seq(-1, 1, length.out = 102)

# Save correlation heatmaps 
save_cor_heatmap <- function(mat, file) {
  mat2 <- na.omit(order_by_mean_minus_sd(as.matrix(mat)))
  pheatmap(mat2, color = my_colors, breaks = my_breaks,
           cluster_rows = FALSE, cluster_cols = FALSE)
  ggsave(file, width = 10, height = 15, dpi = 300)
}

save_cor_heatmap(cor_mat_A, "/vol/output/2025final/stem_signature_correlation_union.pdf") #figure 2b
save_cor_heatmap(cor_mat_B, "/vol/output/2025final/stem_signature_correlation_other_species_union.pdf") #Extended Data Fig. 3b

# Generate signature score scatter plots and per-species z-score heatmaps 
walk(names(objs2), function(sp) {
  sc <- signature_scatter(objs2[[sp]], score_col = "stem1", ylim = c(-0.5, 1))
  # figure 2c: per-species scatter plots
  ggsave(sc$plot, filename = sprintf("/vol/output/2025final/stem_signature_%s_scatter.pdf", sp),
         width = 5, height = 3, dpi = 300)

  genes_for_heat <- mapped_stem[[sp]]
  if (length(genes_for_heat) > 0) {
    # Preserve mouse-based ordering after mapping
    if (!is.null(orth[[sp]])) {
      df <- orth[[sp]]$df[, c(orth[[sp]]$from_col, orth[[sp]]$to_col)] %>% distinct()
      order_vec <- df %>% filter(.data[[orth[[sp]]$from_col]] %in% mouse_stem) %>%
        mutate(mouse_order = match(.data[[orth[[sp]]$from_col]], mouse_stem)) %>%
        arrange(mouse_order) %>%
        pull(.data[[orth[[sp]]$to_col]]) %>% intersect(genes_for_heat)
    } else {
      order_vec <- intersect(mouse_stem, genes_for_heat)
    }

    if (length(order_vec) > 0) {
      pheatmap_obj <- plot_heatmap_simple(objs2[[sp]], order_vec, sc$cell_order)
      # figure 2c: per-species heatmaps
      ggsave(filename = sprintf("/vol/output/2025final/stem_signature_%s_heatmap.pdf", sp),
             plot = pheatmap_obj, width = 5, height = 3, dpi = 300)
    }
  }
})

###################### revival stem cell signature gene score ###########################
mouse_rev <- unique(rev.sig)
# Map mouse revival signature to each species
mapped_rev <- imap(objs, ~ map_genes_to_species(mouse_rev, species = .y, orth_info = orth[[.y]], obj = .x))
# AddModuleScore for each species 
objs3 <- imap(objs, ~ add_signature_score(.x, mapped_rev[[.y]], name = "rev", assay = "RNA"))
# UMAP FeaturePlot with title including the number of signature genes
walk(names(objs3), function(sp) {
  n_genes <- length(mapped_rev[[sp]])
  p <- plot_featureplot_v4(objs3[[sp]], "rev1", 0.8) +
    labs(title = sprintf("Revival signature in %s (#: %d)", sp, n_genes)) +
    theme(plot.title = element_text(size = 20, face = "bold"))
  # Saved as per original: per-species UMAP PDFs
  ggsave(p, filename = sprintf("/vol/output/2025final/revival_signature_%s_umap.pdf", sp),
         width = 5, height = 3, dpi = 300) #Extended Data Fig. 3d
})

# Calculate correlation (score vs. expression) per species
cors_rev <- imap(objs3, ~ cor_vs_score(.x, score_col = "rev1", genes = mapped_rev[[.y]], method = "spearman"))

# Combine species into two panels (as in original union/other species split)
cor_mat <- bind_cor_tables(cors_rev[c("mms","mh","rat","spiny","ferret","bovine","chick","hmouse","hamster","human","marmoset","beagle","pig","goat","hypsugo","steppe","sakhalin")], mouse_ref_genes = mouse_rev)


# Save correlation heatmaps 
save_cor_heatmap_rev <- function(mat, file,
                                 colors = colorRampPalette(c("#213156","#5790BF","white","#C96652","#5F1120"))(101),
                                 breaks = seq(-1, 1, length.out = 102)) {
  # Accept matrix or data.frame
  mat2 <- order_rows_mean_minus_sd_with_NA(mat)
  # pheatmap은 matrix 필요
  mat2_num <- as.data.frame(mat2)
  num_cols <- vapply(mat2_num, is.numeric, logical(1))
  mat2_num <- as.matrix(mat2_num[, num_cols, drop = FALSE])

  # all-NA rows dropped for plotting safety (선택)
  keep <- rowSums(!is.na(mat2_num)) > 0
  mat2_num <- mat2_num[keep, , drop = FALSE]

  pheatmap::pheatmap(mat2_num,
                     color = colors, breaks = breaks,
                     cluster_rows = FALSE, cluster_cols = FALSE)
  ggsave(file, width = 10, height = 15, dpi = 300)
}

save_cor_heatmap_rev(cor_mat, "/vol/output/2025final/rev_signature_correlation_union.pdf") #Extended Data Fig. 3c

###################### secretory signature gene score ###########################
objs_scored <- imap(objs, function(obj, sp) {
  # map each signature to species
  genes_stem <- map_genes_to_species(unique(stem.sig), species = sp, orth_info = orth[[sp]], obj = obj)
  obj <- add_signature_score(obj, genes_stem, name = "stem", assay = "RNA")         # -> stem1

  genes_rev  <- map_genes_to_species(unique(rev.sig),  species = sp, orth_info = orth[[sp]], obj = obj)
  obj <- add_signature_score(obj, genes_rev,  name = "rev",  assay = "RNA")          # -> rev1

  genes_sec  <- map_genes_to_species(unique(sec.sig),  species = sp, orth_info = orth[[sp]], obj = obj)
  obj <- add_signature_score(obj, genes_sec,  name = "secretory", assay = "RNA")     # -> secretory1

  obj
})

options(repr.plot.width = 20, repr.plot.height = 20)

target_species <- c("mms","rat","spiny","ferret","bovine","chick")

plots <- map(target_species, function(sp) {
  obj <- objs_scored[[sp]]
  df_scaled <- prep_radar_data(obj,
                               score_cols  = c("stem1","secretory1","rev1"),
                               cluster_col = "celltype",
                               legend_order = names(my.cols))
  plot_radar(df_scaled, title = sp, colors = my.cols)
})

# NULL 제거
plots <- discard(plots, is.null)

final_plot <- wrap_plots(plots, nrow = 3, ncol = 2)
ggsave("/vol/output/2025final/radar_plot_all_species.pdf",
       final_plot, width = 20, height = 20) #figure 2a

