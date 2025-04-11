#### For generating correlation plots shown in Figures 3 & 4 using scRNA-seq dataset GSE70657

## load working directory and R packages
library(scran)
library(Seurat)
library(dplyr)
library(umap)
library(cowplot)
library(ggplot2)
library(ggridges)
library(ggpubr)
library(GSEABase)
library(GGally)
library(tidyr)
library(AUCell)
library(reshape2)
library(ggcorrplot)
library(escape)
library(dittoSeq)
library(RColorBrewer)
library(grDevices)
load("GSE70657_Selp_HSC.RData")

##################################################
## data matrix was downloaded from: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE70657; GSE70657
## load data matrix and create seurat object
rawmt = as.matrix(read.table("./GSE70657_Grover.A_et.al_RefSeq.Read.Count.txt",header=TRUE, row.names=1))

## creates a Seurat object (sce.rawmt) using the CreateSeuratObject() function. 
sce.rawmt = CreateSeuratObject(counts = rawmt, project = "aHSC", min.cells = 5, min.features =200)

## performs normalization of the gene expression data in the Seurat object (sce.rawmt) using the NormalizeData() function.
sce.rawmt <- NormalizeData(sce.rawmt)

## Adding metadata
sce.rawmt@meta.data$group = ifelse(grepl("old",rownames(sce.rawmt@meta.data)), "Aged", "Young")
sce.rawmt@meta.data$group = factor(sce.rawmt@meta.data$group, levels=c("Young", "Aged"))
sce.rawmt@meta.data$popul = "HSC"

## clustering, and generate UMAP
sce.rawmt = FindVariableFeatures(sce.rawmt, selection.method = "vst", nfeatures = 5000)

## Scaling data
sce.rawmt = ScaleData(sce.rawmt, verbose = FALSE)

## Running Uniform Manifold Approximation and Projection (UMAP)
sce.rawmt = RunUMAP(sce.rawmt, reduction = "pca", dims = 1:20)

## Finding Neighbors and Clusters
sce.rawmt = FindNeighbors(sce.rawmt, reduction = "pca", dims = 1:20)
sce.rawmt = FindClusters(sce.rawmt, resolution = 0.5)
sce.rawmt = FindClusters(sce.rawmt, resolution = 0.8)

## ssGSEA with 'escape' package to calcualte enrichment socore (ES) of gene sets
geneSet_Selphi <- getGmt("./Selp_related_complete.genesets.mouse.gmt")
ES <- enrichIt(obj = sce.rawmt, gene.sets = geneSet_Selphi, groups =2000, cores = 6)
sce.rawmt <- AddMetaData(sce.rawmt, ES)

################################################
comb_meta = sce.rawmt@meta.data
#################################################
#### For Figures 3h, i, j, k
p1 = ggplot(comb_meta, aes(x = sSvendsen_HSCAging_gt4UP, y = Selp_High)) + geom_bin2d() + paletteer::scale_fill_paletteer_c("grDevices::Purples 3") + geom_smooth(method=lm) + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none", text = element_text(size = 15), axis.text=element_text(colour="black"), aspect.ratio = 1) + xlab("Selphi_up") + ylab("HSC_Aging_up") + coord_fixed()
p2 = ggplot(comb_meta, aes(x = Myeloid_high, y = Selp_High)) + geom_bin2d() + paletteer::scale_fill_paletteer_c("grDevices::Purples 3") + geom_smooth(method=lm) + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none", text = element_text(size = 15), axis.text=element_text(colour="black"), aspect.ratio = 1) + xlab("Selphi_up") + ylab("Myeloid-bias") + coord_fixed()
p3 = ggplot(comb_meta, aes(x = IL1B_Act, y = Selp_High)) + geom_bin2d() + paletteer::scale_fill_paletteer_c("grDevices::Purples 3") + geom_smooth(method=lm) + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none", text = element_text(size = 15), axis.text=element_text(colour="black"), aspect.ratio = 1) + xlab("Selphi_up") + ylab("IL-1b_Activation") + coord_fixed()
p4 = ggplot(comb_meta, aes(x = TNF_Act, y = Selp_High)) + geom_bin2d() + paletteer::scale_fill_paletteer_c("grDevices::Purples 3") + geom_smooth(method=lm) + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none", text = element_text(size = 15), axis.text=element_text(colour="black"), aspect.ratio = 1) + xlab("Selphi_up") + ylab("TNF_Activation") + coord_fixed()


pdf("GSE70657.DensityPlot.GSEA.Fig3.pdf", height=5, width=10)
ggarrange(p1, p2, p3, p4, ncol = 4, nrow = 1, align = "v")
dev.off()

cor.test(comb_meta$sSvendsen_HSCAging_gt4UP, comb_meta$Selp_High)
# p-value = 1.0e-08; R = 0.47

cor.test(comb_meta$Myeloid_high, comb_meta$Selp_High)
# p-value = 1.6e-09; R = 0.49

cor.test(comb_meta$Myeloid_high, comb_meta$Selp_High)
# p-value = 2.3e-07; R = 0.43

cor.test(comb_meta$TNF_Act, comb_meta$Selp_High)
# p-value = 1.3e-09; R = 0.49


#### For Figures 4h, i, m
p5 = ggplot(comb_meta, aes(x = IL1B_Act, y = AgedHSC_TXP_Dn)) + geom_bin2d() + paletteer::scale_fill_paletteer_c("grDevices::Purples 3") + geom_smooth(method=lm) + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none", text = element_text(size = 15), axis.text=element_text(colour="black"), aspect.ratio = 1) + xlab("IL-1b_Activation") + ylab("Post-TXP_down") + coord_fixed()
p6 = ggplot(comb_meta, aes(x = TNF_Act, y = AgedHSC_TXP_Dn)) + geom_bin2d() + paletteer::scale_fill_paletteer_c("grDevices::Purples 3") + geom_smooth(method=lm) + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none", text = element_text(size = 15), axis.text=element_text(colour="black"), aspect.ratio = 1) + xlab("TNF_Activation") + ylab("Post-TXP_down") + coord_fixed()
p7 = ggplot(comb_meta, aes(x = sSvendsen_HSCAging_gt4UP, y = AgedHSC_TXP_Dn)) + geom_bin2d() + paletteer::scale_fill_paletteer_c("grDevices::Purples 3") + geom_smooth(method=lm) + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none", text = element_text(size = 15), axis.text=element_text(colour="black"), aspect.ratio = 1) + xlab("HSC_Aging_up") + ylab("Post-TXP_down") + coord_fixed()

pdf("GSE70657.DensityPlot.GSEA.Fig4.pdf", height=5, width=10)
ggarrange(p5, p6, p7, ncol = 3, nrow = 1, align = "v")
dev.off()

cor.test(comb_meta$IL1B_Act, comb_meta$AgedHSC_TXP_Dn)
# p-value = 2.2e-16; R = 0.69

cor.test(comb_meta$TNF_Act, comb_meta$AgedHSC_TXP_Dn)
# p-value = 2.2e-16; R = 0.77

cor.test(comb_meta$sSvendsen_HSCAging_gt4UP, comb_meta$AgedHSC_TXP_Dn)
# p-value = 2.2e-16; R = 0.88












