library(pheatmap)
library(dplyr)
library(readr)
library(tibble)
library(edgeR)
library(tximport)
library(writexl)
library(clusterProfiler)

root_path = "./"

# result path
result_path = paste0(root_path,"Fig_4_heatmap/")
if (!file.exists(result_path)) {
  dir.create(result_path)
}

SPINY_YAP = read_csv("./Org_Zoo_Bulk_RNA_seq_analysis/SMP_SMC_TPM.csv")
MOUSE_YAP = read_csv("./Org_Zoo_Bulk_RNA_seq_analysis/MP_MC_TPM.csv")
SPINY_Ble = read_csv("./Org_Zoo_Bulk_RNA_seq_analysis/SP_B_SP_C_TPM.csv")
MOUSE_Ble = read_csv("./Org_Zoo_Bulk_RNA_seq_analysis/MO_B_MO_C_TPM.csv")

# Remove unnecessary columns
drop_cols <- c("...1", "gene_id")
MOUSE_YAP <- MOUSE_YAP %>% select(-any_of(drop_cols))
SPINY_YAP <- SPINY_YAP %>% select(-any_of(drop_cols))
MOUSE_Ble <- MOUSE_Ble %>% select(-any_of(drop_cols))
SPINY_Ble <- SPINY_Ble %>% select(-any_of(drop_cols))

# Load gene sets
geneset_gmt = read.gmt("./Bulk_RNA_seq_geneset.gmt")
P53_target_visualization_genes <- geneset_gmt %>%
  filter(term == "P53_target_visualization") %>%
  pull(gene)
Yap_target_visualization_genes <- geneset_gmt %>%
  filter(term == "Yap_target_visualization") %>%
  pull(gene)


## <Bleomycin case - Fig4k>
# Prepare matrix for heatmap
MOUSE_Ble = MOUSE_Ble[MOUSE_Ble$gene_name %in% P53_target_visualization_genes, ]
SPINY_Ble = SPINY_Ble[SPINY_Ble$gene_name %in% P53_target_visualization_genes, ]
TPM_mat_Ble <- MOUSE_Ble %>%
  inner_join(SPINY_Ble, by = "gene_name") %>%
  as.data.frame() %>%
  column_to_rownames("gene_name") %>%
  as.matrix()

# Generate heatmap
log2_tpm_Ble = log2(TPM_mat_Ble + 1)
zscore_tpm_Ble <- t(scale(t(log2_tpm_Ble)))

plot_tpm_Ble = zscore_tpm_Ble[,c("MO_C1","MO_C2","MO_C3","MO_B1","MO_B2","MO_B3","SP_C1","SP_C2","SP_C3","SP_B1","SP_B2","SP_B3")]
plot_tpm_Ble = plot_tpm_Ble[P53_target_visualization_genes,]
outfile = paste0(result_path,"Fig4k_heatmap_Bleomycin.pdf")
heatmap_plot = pheatmap(plot_tpm_Ble,
         scale = "none",   
         filename = outfile,
         show_rownames = TRUE,
         border_color = NA,
         cluster_rows = FALSE,
         cluster_cols = FALSE
         )

## <PY-60 case - Fig4f>
# Prepare matrix for heatmap
MOUSE_YAP = MOUSE_YAP[MOUSE_YAP$gene_name %in% Yap_target_visualization_genes, ]
SPINY_YAP = SPINY_YAP[SPINY_YAP$gene_name %in% Yap_target_visualization_genes, ]
TPM_mat_YAP <- MOUSE_YAP %>%
  inner_join(SPINY_YAP, by = "gene_name") %>%
  as.data.frame() %>%
  column_to_rownames("gene_name") %>%          
  as.matrix()

log2_tpm_YAP = log2(TPM_mat_YAP + 1)
zscore_tpm_YAP <- t(scale(t(log2_tpm_YAP)))

# Generate heatmap
plot_tpm_YAP = zscore_tpm_YAP[,c("MC1", "MC2", "MC3","MP1", "MP2", "MP3", "SMC1", "SMC2", "SMC3", "SMP1", "SMP2", "SMP3")]
plot_tpm_YAP = plot_tpm_YAP[Yap_target_visualization_genes,]
outfile = paste0(result_path,"Fig4f_heatmap_Yap.pdf")
heatmap_plot = pheatmap(plot_tpm_YAP,
         scale = "none",   
         filename = outfile,
         show_rownames = TRUE,
         border_color = NA,
         cluster_rows = FALSE,
         cluster_cols = FALSE
         )