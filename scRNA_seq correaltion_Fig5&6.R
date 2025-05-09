#### For generating correlation plots shown in Figures 5&6 using scRNA-seq dataset GSE232022

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

load("Selp-CITEseq-IC.RData")

## load expression matrix generated by Cellranger

ct_dir = './Control/sample_filtered_feature_bc_matrix'
list.files(ct_dir)
ct.data <- Read10X(data.dir = ct_dir)

dfo_dir = './DFO/sample_filtered_feature_bc_matrix'
list.files(dfo_dir)
dfo.data <- Read10X(data.dir = dfo_dir)


## Analysis of CITE-seq data was performed according to the seurat manuel: https://satijalab.org/seurat/articles/multimodal_vignette
## creates a Seurat object
ct = CreateSeuratObject(counts = ct.data$`Gene Expression`)
ct[['Protein']] = CreateAssayObject(counts = ct.data$`Antibody Capture`)
ct[['Hashtag']] = CreateAssayObject(counts = ct.data$`Multiplexing Capture`)

dfo = CreateSeuratObject(counts = dfo.data$`Gene Expression`)
dfo[['Protein']] = CreateAssayObject(counts = dfo.data$`Antibody Capture`)
dfo[['Hashtag']] = CreateAssayObject(counts = dfo.data$`Multiplexing Capture`)

ct[["percent.mt"]] = PercentageFeatureSet(ct, pattern = "^mt-")
dfo[["percent.mt"]] = PercentageFeatureSet(dfo, pattern = "^mt-")

rm(ct.data, dfo.data)


## merge the objects of DFO and Control
ct$treatment = "Control"
dfo$treatment = "DFO"

## merge the seurat objects of NT and DFO, and filter out cells with low sequencing quality, dead cells, or potential doublets
combined <- merge(ct, y = dfo, add.cell.ids = c("ct", "dfo"), project = "IC_scHSC")
combined <- subset(combined, subset = nFeature_RNA > 500 & nFeature_RNA < 8000 & percent.mt < 15 & nCount_RNA < 100000 & nCount_RNA > 2000)

## normalization for RNA and protein data are performed seperately
combined <- NormalizeData(combined)
combined <- NormalizeData(combined, assay = "Protein", normalization.method = "CLR")



## Identify HSC based on the protein expression (presented by levels of antibody-derived tags (ADT)) of CD150 and CD48
## Check the expression of CD150 and CD48, with pseudo-gating
DefaultAssay(combined) <- "Protein"
library(scales)

proteins = c("Ms.CD150","Ms.CD48")
pro_slam = GetAssayData(object = combined, slot = 'data')[proteins, ]
pro_slam = t(pro_slam)

## add the ADT level of SLAM markers, CD150 and CD48, to meta.data
combined@meta.data$ADT_CD150 = GetAssayData(object = combined, assay = 'Protein')['Ms.CD150', ]
combined@meta.data$ADT_CD48 = GetAssayData(object = combined, assay = 'Protein')['Ms.CD48', ]

## HSC is defined with SLAM gating, CD150+CD48-
combined@meta.data$HSPC = ifelse((combined@meta.data$ADT_CD150>1 & combined@meta.data$ADT_CD48<1),"HSC",ifelse((combined@meta.data$ADT_CD150<1 & combined@meta.data$ADT_CD48<1),"MPP",ifelse((combined@meta.data$ADT_CD150>1 & combined@meta.data$ADT_CD48>1),"HPC","Prog")))
combined@meta.data$HSPC =  as.factor(combined@meta.data$HSPC)
table(combined@meta.data$HSPC)

## subset HSC cells
combined_hsc <- subset(combined, subset = HSPC == "HSC")

##################################################
## Run the standard workflow for visualization and clustering
## In this study, the library prepration and sequencing of Control or DFO samples were performed in the same batch
## If not, integration analysis is needed with the merged object before clustering
## tutorial: https://satijalab.org/seurat/articles/integration_introduction.html

DefaultAssay(combined_hsc) <- "RNA"

combined_hsc <- FindVariableFeatures(combined_hsc)
combined_hsc <- ScaleData(combined_hsc)
combined_hsc <- RunPCA(combined_hsc, verbose = FALSE)

combined_hsc <- FindNeighbors(combined_hsc, dims = 1:25)
combined_hsc <- RunUMAP(combined_hsc, dims = 1:25)
combined_hsc <- RunTSNE(combined_hsc, dims = 1:25, method = "FIt-SNE")
combined_hsc <- FindClusters(combined_hsc, resolution = 0.8, verbose = FALSE)

ident_colours = c('#F8766D', '#D39200', '#93AA00', '#00BA38', '#00C19F', '#00B9E3', '#619CFF', '#DB72FB', '#FF61C3')

## UMAP dim plot, splitted by treatment conditions, or both conditions together
pdf("Seurat.HSConly.UMAP.combine.res08.pdf", height=4, width=4)
DimPlot(combined_hsc, reduction = "umap", label = TRUE, label.size = 5, cols = ident_colours)
dev.off()

pdf("Seurat.HSConly.UMAP.res08.pdf", height=4, width=6)
DimPlot(combined_hsc, reduction = "umap", split.by = "treatment", label = TRUE, label.size = 5)
dev.off()

## check the percentage of mt-DNA in each cluster
pdf("20230413_Percent-mt.pdf", height=4, width=6)
VlnPlot(combined_hsc, features="percent.mt", y.max=20, ncol = 1,cols = col)
dev.off()

## we noticed cluster 8 has much higher mt-DNA compared to other clusters, excluded from downstream analysis
combined_hsc <- subset(combined_hsc, idents = c(0, 1, 2, 3, 4, 5, 6,7))

###################################################
## ssGSEA with 'escape' package to calcualte enrichment socore (ES) of gene sets
## GSEA should be performed with RNA data, not protein or ADT.
DefaultAssay(combined_hsc) <- "RNA"

geneSet_Selphi <- getGmt("./Selp_related_complete.genesets.mouse.gmt")
ES <- enrichIt(obj = combined_hsc, gene.sets = geneSet_Selphi, groups =2000, cores = 10)
combined_hsc <- AddMetaData(combined_hsc, ES)


################################################
comb_meta = combined_hsc@meta.data
#################################################
#### For Figures 5f,g
p1 = ggplot(comb_meta, aes(x = PSGL1_down, y = AgedHSC_TXP_Dn)) + geom_bin2d() + paletteer::scale_fill_paletteer_c("grDevices::Purples 3") + geom_smooth(method=lm) + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none", text = element_text(size = 15), axis.text=element_text(colour="black"), aspect.ratio = 1) + xlab("PSGL1_down") + ylab("Post-TXP_down") + coord_fixed()
p2 = ggplot(comb_meta, aes(x = PSGL1_down, y = sSvendsen_HSCAging_gt4UP)) + geom_bin2d() + paletteer::scale_fill_paletteer_c("grDevices::Purples 3") + geom_smooth(method=lm) + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none", text = element_text(size = 15), axis.text=element_text(colour="black"), aspect.ratio = 1) + xlab("PSGL1_down") + ylab("HSC_Aging_up") + coord_fixed()

pdf("GSE232022.DensityPlot.GSEA.Fig5.pdf", height=4, width=6)
ggarrange(p3, p4, p5, ncol = 3, nrow = 1, align = "v")
dev.off()

cor.test(comb_meta$PSGL1_down, comb_meta$AgedHSC_TXP_Dn)
# p-value < 2.2e-16; R = 0.58

cor.test(comb_meta$PSGL1_down, comb_meta$sSvendsen_HSCAging_gt4UP)
# p-value < 2.2e-16; R = 0.26



#### For Figures 6e, f, g
p3 = ggplot(comb_meta, aes(x = SelpKO_Up, y = sSvendsen_HSCAging_gt4UP)) + geom_bin2d() + paletteer::scale_fill_paletteer_c("viridis::plasma") + geom_smooth(method=lm) + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none", text = element_text(size = 15), axis.text=element_text(colour="black"), aspect.ratio = 1) + xlab("SelpKO_up") + ylab("HSC_Aging_up") + scale_color_manual(values=ident_colours) + coord_fixed()
p4 = ggplot(comb_meta, aes(x = SelpKO_Up, y = IL6_Act)) + geom_bin2d() + paletteer::scale_fill_paletteer_c("viridis::plasma") + geom_smooth(method=lm) + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none", text = element_text(size = 15), axis.text=element_text(colour="black"), aspect.ratio = 1) + xlab("SelpKO_up") + ylab("IL6_Activation") + scale_color_manual(values=ident_colours) + coord_fixed()
p5 = ggplot(comb_meta, aes(x = SelpKO_Up, y = IL1B_Act)) + geom_bin2d() + paletteer::scale_fill_paletteer_c("viridis::plasma") + geom_smooth(method=lm) + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none", text = element_text(size = 15), axis.text=element_text(colour="black"), aspect.ratio = 1) + xlab("SelpKO_up") + ylab("IL-1b_Activation") + scale_color_manual(values=ident_colours) + coord_fixed()

pdf("GSE232022.DensityPlot.GSEA.Fig6.pdf", height=4, width=6)
ggarrange(p3, p4, p5, ncol = 3, nrow = 1, align = "v")
dev.off()

cor.test(comb_meta$SelpKO_Up, comb_meta$sSvendsen_HSCAging_gt4UP)
# p-value < 2.2e-16; R = 0.62 

cor.test(comb_meta$SelpKO_Up, comb_meta$IL6_Act)
# p-value < 2.2e-16; R = 0.52 

cor.test(comb_meta$SelpKO_Up, comb_meta$IL1B_Act)
# p-value < 2.2e-16; R = 0.40













