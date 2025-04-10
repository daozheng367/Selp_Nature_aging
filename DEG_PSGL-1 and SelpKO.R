
library(DESeq2)
library(pheatmap)
library(RColorBrewer)
library(stringr)
tab <- read.csv('Daozheng_BulkRNAseq2_FeatureCount_Raw.csv')
cts <- as.matrix(tab[,-1])
rownames(cts) <- tab$X
SampleName <- colnames(cts)
Condition <- c(rep('KO', 3), rep('KO_LG', 3), rep('KO_VH', 3),
               rep('WT', 3), rep('WT_LG', 3), rep('WT_VH', 3))
coldata <- data.frame(SampleName = SampleName,
                      condition = Condition)
coldata$condition <- factor(coldata$condition)
# 
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ condition)
dds
# DEG analysis
dds <- dds[rowSums(counts(dds))>10, ]
dds <- DESeq(dds)
rld <- rlog(dds, blind = FALSE)
sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix, clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists, cluster_cols=F, col=colors)
plotPCA(rld)

pca_df <- plotPCA(rld, intgroup='condition', returnData=T)
library(ggplot2)
names <- rownames(sampleDistMatrix)
ggplot(pca_df, aes(x=PC1, y=PC2, color=condition)) +
  geom_point() + 
  geom_text(aes(label = names), hjust = 1, vjust = 0, nudge_x = -0.1, nudge_y = 0.1) 

### Comparison1:  WT vs KO ===========================================================
head(tab)
tab1 <- tab[ , c(2:4, 11:13)]
head(tab1)
cts1 <- as.matrix(tab1)
rownames(cts1) <- tab$X
head(cts1)
SampleName1 <- colnames(cts1)
Condition1 <- c(rep('KO', 3), rep('WT', 3))
coldata1 <- data.frame(SampleName = SampleName1,
                      condition = Condition1)
coldata1$condition <- factor(coldata1$condition)
# 
dds1 <- DESeqDataSetFromMatrix(countData = cts1,
                              colData = coldata1,
                              design = ~ condition)
dds1
# DEG analysis
dds1 <- dds1[rowSums(counts(dds1))>10, ]
dds1$condition <- factor(dds1$condition, levels = c('WT', 'KO'))
dds1 <- DESeq(dds1)
rld1 <- rlog(dds1, blind = FALSE)
sampleDists1 <- dist(t(assay(rld1)))
sampleDistMatrix1 <- as.matrix(sampleDists1)
colnames(sampleDistMatrix1) <- NULL
pheatmap(sampleDistMatrix1, clustering_distance_rows = sampleDists1,
         clustering_distance_cols = sampleDists1, cluster_cols=F, col=colors)
plotPCA(rld1)
### read counts
Count_dds1 <- counts(dds1, normalized=T)
res1 <- results(dds1)
res1
# output spreadsheet
dat1 <- cbind(res1, Count_dds1)
write.csv(dat1, 'DEGs_KO_WT_RNAseq.csv')

summary(res1)
res1[order(res1$padj),] 
resLFC1 <- lfcShrink(dds1, coef='condition_KO_vs_WT', type='apeglm')
plotMA(resLFC1, main="KO vs. WT", ylim=c(-6,16))

### Comparison1:  WT  LG vs VH ===========================================================
head(tab)
tab2 <- tab[ , c(14:19)]
tab2 <- tab[ , c(15:19)]  # take WT_LG15 out
head(tab2)
cts2 <- as.matrix(tab2)
rownames(cts2) <- tab$X
head(cts2)
SampleName2 <- colnames(cts2)
Condition2 <- c(rep('WT_LG', 2), rep('WT_VH', 3))
coldata2 <- data.frame(SampleName = SampleName2,
                       condition = Condition2)
coldata2$condition <- factor(coldata2$condition)
# 
dds2 <- DESeqDataSetFromMatrix(countData = cts2,
                               colData = coldata2,
                               design = ~ condition)
dds2
# DEG analysis
dds2 <- dds2[rowSums(counts(dds2))>10, ]
dds2$condition <- factor(dds2$condition, levels = c('WT_VH', 'WT_LG'))
dds2 <- DESeq(dds2)
rld2 <- rlog(dds2, blind = FALSE)
sampleDists2 <- dist(t(assay(rld2)))
sampleDistMatrix2 <- as.matrix(sampleDists2)
colnames(sampleDistMatrix2) <- NULL
pheatmap(sampleDistMatrix2, clustering_distance_rows = sampleDists2,
         clustering_distance_cols = sampleDists2, cluster_cols=F, col=colors)
plotPCA(rld2)
### read counts
Count_dds2 <- counts(dds2, normalized=T)
res2 <- results(dds2)
res2
# output spreadsheet
dat2 <- cbind(res2, Count_dds2)
write.csv(dat2, 'DEGs_WTLG_WTVH_RNAseq.csv')

summary(res2)
res2[order(res2$padj),] 
resLFC2 <- lfcShrink(dds2, coef='condition_WT_LG_vs_WT_VH', type='apeglm')
plotMA(resLFC2, main="WT_LG vs. WT_VH", ylim=c(-16,8))

### Comparison1:  KO  LG vs VH ===========================================================
head(tab)
tab3 <- tab[ , c(5:10)]
tab3 <- tab[ , c(5,7:9)] # take KO_LG18 and KO_VH5
head(tab3)
cts3 <- as.matrix(tab3)
rownames(cts3) <- tab$X
head(cts3)
SampleName3 <- colnames(cts3)
Condition3 <- c(rep('KO_LG', 2), rep('KO_VH', 2))
coldata3 <- data.frame(SampleName = SampleName3,
                       condition = Condition3)
coldata3$condition <- factor(coldata3$condition)
# 
dds3 <- DESeqDataSetFromMatrix(countData = cts3,
                               colData = coldata3,
                               design = ~ condition)
dds3
# DEG analysis
dds3 <- dds3[rowSums(counts(dds3))>10, ]
dds3$condition <- factor(dds3$condition, levels = c('KO_VH', 'KO_LG'))
dds3 <- DESeq(dds3)
rld3 <- rlog(dds3, blind = FALSE)
sampleDists3 <- dist(t(assay(rld3)))
sampleDistMatrix3 <- as.matrix(sampleDists3)
colnames(sampleDistMatrix3) <- NULL
pheatmap(sampleDistMatrix3, clustering_distance_rows = sampleDists3,
         clustering_distance_cols = sampleDists3, cluster_cols=F, col=colors)
plotPCA(rld3)

pca_df <- plotPCA(rld3, intgroup='condition', returnData=T)
library(ggplot2)
names <- rownames(sampleDistMatrix3)
ggplot(pca_df, aes(x=PC1, y=PC2, color=condition)) +
  geom_point() + 
  geom_text(aes(label = names), hjust = 1, vjust = 0, nudge_x = -0.1, nudge_y = 0.1) 

### read counts
Count_dds3 <- counts(dds3, normalized=T)
res3 <- results(dds3)
res3
# output spreadsheet
dat3 <- cbind(res3, Count_dds3)
write.csv(dat3, 'DEGs_KOLG_KOVH_RNAseq.csv')

summary(res3)
res3[order(res3$padj),] 
resLFC3 <- lfcShrink(dds3, coef='condition_KO_LG_vs_KO_VH', type='apeglm')
plotMA(resLFC3, main="KO_LG vs. WT_VH", ylim=c(-5,15))






