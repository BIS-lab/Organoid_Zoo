library(pheatmap)
library(dplyr)
library(readr)
library(tibble)
library(edgeR)
library(tximport)
library(writexl)
library(clusterProfiler)

root_path = "./"

save_edgeR_DEG = function(path, data, Comparison_name){
    result_data = data %>% filter(!is.na(FDR))
    dge_description_df <- data.frame(
    Column = c("logCPM", "logFC", "F", "PValue", "FDR", "gene_id", "gene_name"),
    Meaning = c(
        "Log Normalized average expression value across all samples in two groups",
        "Log2-transformed fold change between two groups",
        "F-statistic from quasi-likelihood test from EdgeR",
        "Raw enrichment p-value with quasi-likelihood F-tests",
        "Adjusted p-value (false discovery rate (FDR))",
        "Ensembl gene id",
        "Gene Symbol"
    )
    )

    # write DGE data
    result_path = paste0(path,"0.DGE/")
    if (!file.exists(result_path)) {
      dir.create(result_path)
    }
    out_list <- list(
        "Result" = result_data,
        "Column_Description" = dge_description_df
    )
    write_xlsx(out_list, paste0(result_path,"DEG_edgeR_",Comparison_name,".xlsx"))
}
 

project_settings <- list(

  list(
    Project_name = "Org_Zoo_Bulk_RNA_seq_analysis/",
    comb_set = list(c("MP", "MC")),
    tx2gene_path = "./0_result_mouse/star_salmon/tx2gene.tsv",
    metadata_path = "./Metadata/metadata_mouse_YAP.csv"
  ),

  list(
    Project_name = "Org_Zoo_Bulk_RNA_seq_analysis/",
    comb_set = list(c("MO_B", "MO_C")),
    tx2gene_path = "./0_result_mouse/star_salmon/tx2gene.tsv",
    metadata_path = "./Metadata/metadata_mouse_Ble.csv"
  ),

  list(
    Project_name = "Org_Zoo_Bulk_RNA_seq_analysis/",
    comb_set = list(c("SMP", "SMC")),
    tx2gene_path = "./0_result_spiny/star_salmon/tx2gene.tsv",
    metadata_path = "./Metadata/metadata_spiny_YAP.csv"
  ),

  list(
    Project_name = "Org_Zoo_Bulk_RNA_seq_analysis/",
    comb_set = list(c("SP_B", "SP_C")),
    tx2gene_path = "./0_result_spiny/star_salmon/tx2gene.tsv",
    metadata_path = "./Metadata/metadata_spiny_Ble.csv"
  )
)

# make result file
DGE_result = list()
      
for (setting in project_settings) {
  Project_name <- setting$Project_name
  comb_set <- setting$comb_set
  tx2gene_path <- setting$tx2gene_path
  metadata_path <- setting$metadata_path
  
  # Create project directory
  result_path = paste0(root_path,Project_name)
  if (!file.exists(result_path)) {
    dir.create(result_path)
  }

  # load count_data & col_data & description
  tx2gene_data = as.data.frame(read.csv(tx2gene_path ,sep="\t"))
  Col_data = as.data.frame(read.csv(metadata_path))
  files <- setNames(Col_data$file, Col_data$sample)
  stopifnot(all(file.exists(files)))

  txi.salmon <- tximport(files, type = "salmon", tx2gene = tx2gene_data)
  Col_data <- Col_data[match(colnames(txi.salmon$counts), Col_data$sample), ]


  # gene_id to gene_name template
  tx2gene_unique <- tx2gene_data %>%
    dplyr::select(gene_id, gene_name) %>%
    distinct(gene_id, .keep_all = TRUE)

  ## Save TPM matrix
  tpm_matrix <- txi.salmon$abundance
  tpm_data <- tpm_matrix %>%
    as.data.frame() %>%
    rownames_to_column(var = "gene_id") %>%
    left_join(tx2gene_unique, by = "gene_id") %>%
    relocate(gene_id, .before = 1) %>%
    relocate(gene_name, .after = gene_id)
  tpm_final<- as.matrix(tpm_data)
  write.csv(tpm_final, file = paste0(result_path,comb_set[[1]][1],"_",comb_set[[1]][2],"_TPM.csv"))


  cts <- txi.salmon$counts
  normMat <- txi.salmon$length

  # Obtaining per-observation scaling factors for length, adjusted to avoid
  # changing the magnitude of the counts.
  normMat <- normMat/exp(rowMeans(log(normMat))) # geometric mean of length factor across the samples
  normCts <- cts/normMat

  # Computing effective library sizes from scaled counts, to account for
  # composition biases between samples.
  eff.lib <- calcNormFactors(normCts) * colSums(normCts)

  # Combining effective library sizes with the length factors, and calculating
  # offsets for a log-link GLM.
  normMat <- sweep(normMat, 2, eff.lib, "*")
  normMat <- log(normMat)

  # Creating a DGEList object for use in edgeR.
  y <- DGEList(cts)
  y <- scaleOffset(y, normMat)

  design <-model.matrix(~ 0 + Condition, data = Col_data)
  colnames(design) <- gsub("^Condition", "", colnames(design))

  keep <-filterByExpr(y, design)
  y <- y[keep, ]

  dge_all <- estimateDisp(y, design)
  fit_all <- glmQLFit(dge_all, design)

  for (pair in comb_set) {
    treat <- pair[1]
    con <- pair[2]
    key_name <- paste0(treat, "_vs_", con)

    comparison_set = paste0(treat,"-",con)
    contrast_df = makeContrasts(contrasts = comparison_set,levels=design)


    qlf <- glmQLFTest(fit_all, contrast = contrast_df)
    result_df <- as.data.frame(topTags(qlf, n = Inf, adjust.method = "BH"))
    
    result_df$gene_id <- rownames(result_df)
    result_df <- result_df[, c("gene_id", setdiff(colnames(result_df), "gene_id"))]
    result_df = result_df %>%
      left_join(tx2gene_unique, by = "gene_id")

    DGE_result[[key_name]] <- result_df

  }

  for (Comparison_name in names(DGE_result)){
      save_edgeR_DEG(path = result_path,
                      data = DGE_result[[Comparison_name]],
                      Comparison_name = Comparison_name
                      )
  } 

}