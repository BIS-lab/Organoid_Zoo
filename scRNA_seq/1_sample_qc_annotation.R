## As a representative example, we demonstrate the analysis process for mouse IcM medium organoids, and the analysis of other species also followed the same process.
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(DropletUtils)
library(DoubletFinder)
library(SeuratDisk)
library(biomaRt)

### Cellranger version 7.2.0 was used for the alignment of raw data to the reference genome and generation of count matrix.
### load raw data and create Seurat object
mms.data <- Read10X(data.dir = "/vol/data/raw_Data/mms/outs/raw_feature_bc_matrix/")
mms <- CreateSeuratObject(counts = mms.data, project = "mouse_ms")

##### emptyDrops filtering
temp <- mms[["RNA"]]@counts
out <- emptyDrops(temp)
is.cell <- out$FDR <= 0.001
sum(is.cell, na.rm=TRUE)

mms.counts <- temp[,which(is.cell),drop=FALSE]

mms <- CreateSeuratObject(counts = mms.counts, project = "mouse_ms", min.cells = 3, min.features = 200)
mms[["percent.mt"]] <- PercentageFeatureSet(mms, pattern = "^mt-")

#### QC filtering, For QC information for other species, please refer to the extended data table.
mms <- subset(mms, subset = nFeature_RNA > 1900 & nCount_RNA < 50000 & percent.mt < 15)

s.genes.human <- cc.genes.updated.2019$s.genes
g2m.genes.human <- cc.genes.updated.2019$g2m.genes

# human to mouse ortholog mapping
ortholog_table_filtered <- subset(ortholog_table, hsapiens_homolog_associated_gene_name %in% c(s.genes.human, g2m.genes.human))
ortholog_table_filtered <- ortholog_table_filtered[!duplicated(ortholog_table_filtered$hsapiens_homolog_associated_gene_name), ]

# S phase genes
f.s.genes <- ortholog_table_filtered$external_gene_name[
  ortholog_table_filtered$hsapiens_homolog_associated_gene_name %in% s.genes.human
]
# G2/M phase genes
f.g2m.genes <- ortholog_table_filtered$external_gene_name[
  ortholog_table_filtered$hsapiens_homolog_associated_gene_name %in% g2m.genes.human
]

mms <- CellCycleScoring(mms, s.features = f.s.genes, g2m.features = f.g2m.genes, set.ident = TRUE)

####Normalization & hvg & scaling
mms <- NormalizeData(mms)
mms <- FindVariableFeatures(mms, selection.method = "vst", nfeatures = 2000)

###cell cycling regressing out
mms$CC.Difference <- mms$S.Score - mms$G2M.Score
mms <- ScaleData(mms, vars.to.regress = "CC.Difference", features = rownames(mms) ,verbose = FALSE)

####visualization and clustering
mms <- RunPCA(mms, npcs = 30, verbose = FALSE)
mms <- RunUMAP(mms, reduction = "pca", dims = 1:30)
mms <- FindNeighbors(mms, reduction = "pca", dims = 1:30)
mms <- FindClusters(mms, algorithm = 1, resolution = 0.6)

#### DoubletFinder
sweep.res.list.mms<-paramSweep(mms, PCs = 1:30, sct = FALSE)
sweep.stats.mms<-summarizeSweep(sweep.res.list.mms, GT = FALSE)
bcmvn.mms<-find.pK(sweep.stats.mms)

homotypic.prop <- modelHomotypic(mms@meta.data$nFeature_RNA)
nExp_poi<-round(0.08*nrow(mms@meta.data))

mms <- doubletFinder(mms, PCs = 1:30, pN = 0.25, pK = 0.25, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE) #select pK = 0.25 from bcmvn.mms
mms <- doubletFinder(mms, PCs = 1:30, pN = 0.25, pK = 0.25, nExp = nExp_poi, reuse.pANN = "DF.classifications_0.25_0.25_489", sct = FALSE) #select reuse.pANN from previous step

mms.singlet <- subset(mms, subset = DF.classifications_0.25_0.25_489 == "Singlet")

####Normalization & hvg & scaling
mms.singlet <- NormalizeData(mms.singlet)
mms.singlet <- FindVariableFeatures(mms.singlet, selection.method = "vst", nfeatures = 2000)

###cell cycling regressing out
mms.singlet$CC.Difference <- mms.singlet$S.Score - mms.singlet$G2M.Score
mms.singlet <- ScaleData(mms.singlet, vars.to.regress = "CC.Difference", features = rownames(mms.singlet) ,verbose = FALSE)

####visualization and clustering
mms.singlet <- RunPCA(mms.singlet, npcs = 30, verbose = FALSE)
mms.singlet <- RunUMAP(mms.singlet, reduction = "pca", dims = 1:30)
mms.singlet <- FindNeighbors(mms.singlet, reduction = "pca", dims = 1:30)
mms.singlet <- FindClusters(mms.singlet, resolution = 1.2) #For resolution parameter for other species, please refer to the extended data table.


### Load functions for average and proportion dotplots
source("function_Averageplot.R")

#####combined cell marker
stem <- c("Lgr5","Axin2","Lrig1","Olfm4","Sp5","Ascl2","Slc12a2")
proliferate <- c("Mki67","Cdk4","Mcm5","Mcm6","Pcna","Lig1","Stmn1")
enterocyte <- c("Arg2","Il18","Ccl25","Alpi","Apoa4","Apoa1","Apoc3","Rbp2","Aldob","Slc2a2")
rev <- c("Clu","Anxa1","Anxa3","Ccnd1","Ccnd2")
gp <- c("Lyz1","Gfi1","Spdef","Tff3","Muc2","Muc4","Agr2","Defa24","Atoh1","Neurog3","Dll1","Clca1")
ee <- c("Chga","Chgb")
Tuft <- c("Dclk1","Trpm5")

## if you want to plot specific species markers, you can use the following code
#example for human
# orth <- read.csv("/vol/data/mms_hsapiens_ortholog.csv", stringsAsFactors = FALSE)
# orth <- subset(orth, hsapiens_homolog_orthology_type == "ortholog_one2one")
# map_gene <- function(g) unique(orth$hsapiens_homolog_associated_gene_name[orth$external_gene_name %in% g])
# stem        <- map_gene(stem)
# proliferate <- map_gene(proliferate)
# enterocyte  <- map_gene(enterocyte)
# rev         <- map_gene(rev)
# gp          <- map_gene(gp)
# ee          <- map_gene(ee)
# Tuft        <- map_gene(Tuft)

options(repr.plot.width=12, repr.plot.height=3)

f4 <- average_dotplot(mms, c(stem,proliferate,enterocyte,rev,gp,ee,Tuft),group.by = "celltype", cl.order = c("Stem cell","Proliferating","Absorptive","Revival stem cell","Secretory")) + NoLegend()
f4
ggsave(plot = f4, width = 12, height = 3, dpi = 300, filename = "/vol/output/2025final/mouse_averageplot_celltype_combine_marker.pdf") #Fig2a, Extended Data Fig2b, Extended Data Fig5c dotplot

mms.singlet$celltype <- mms.singlet$seurat_clusters
levels(mms.singlet$celltype)[match('0',levels(mms.singlet$celltype))] <- "Revival stem cell"
levels(mms.singlet$celltype)[match('1',levels(mms.singlet$celltype))] <- "Stem cell"
levels(mms.singlet$celltype)[match('2',levels(mms.singlet$celltype))] <- "Stem cell"
levels(mms.singlet$celltype)[match('3',levels(mms.singlet$celltype))] <- "Proliferating"
levels(mms.singlet$celltype)[match('4',levels(mms.singlet$celltype))] <- "Proliferating"
levels(mms.singlet$celltype)[match('5',levels(mms.singlet$celltype))] <- "Proliferating"
levels(mms.singlet$celltype)[match('6',levels(mms.singlet$celltype))] <- "Stem cell"
levels(mms.singlet$celltype)[match('7',levels(mms.singlet$celltype))] <- "Stem cell"
levels(mms.singlet$celltype)[match('8',levels(mms.singlet$celltype))] <- "Stem cell"
levels(mms.singlet$celltype)[match('9',levels(mms.singlet$celltype))] <- "Proliferating"
levels(mms.singlet$celltype)[match('10',levels(mms.singlet$celltype))] <- "Secretory"
levels(mms.singlet$celltype)[match('11',levels(mms.singlet$celltype))] <- "Secretory"
levels(mms.singlet$celltype)[match('12',levels(mms.singlet$celltype))] <- "Stem cell"
levels(mms.singlet$celltype)[match('13',levels(mms.singlet$celltype))] <- "Absorptive"

my.cols <- c("Stem cell"="#9A4758","Absorptive"="#F0E58D",'Revival stem cell'="#3E5599",'Proliferating'="#578C7F",
'Secretory'="#9D78AE")

options(repr.plot.width=6, repr.plot.height=6)

p1 <- DimPlot(mms, reduction = "umap",group.by="celltype", label = FALSE, repel = TRUE, pt.size = 1, cols=my.cols) + NoLegend()
p1
ggsave(plot = p1, width = 12, height = 12, dpi = 300, filename = "/vol/output/2025final/mouse_averageplot_celltype_umap.pdf") #Fig2a, Extended Data Fig2a, Fig3e UMAP

### Save final singlet object for future analysis
SaveH5Seurat(mms.singlet, filename = "../data/mouse_Ms_singlet.h5Seurat")
