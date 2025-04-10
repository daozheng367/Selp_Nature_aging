
install.packages("dendsort")

library(pheatmap)
library(dendsort)
library(dplyr)


####################################get reads of selected genes
#get full list pre-selected gene list
overlap_gene<-read.csv("20240509_SelpKO_AS_genes.csv",header=T)
head(overlap_gene)
row.names(overlap_gene)
#get normalized reads count
reads<-read.csv("Reads_KO_WT_RNAseq_noKO4_normalized_.csv",header=T)
head(reads)
###make subset of normalized expressions of overlap_gene
reads_eoverlap_gene<-subset(reads, reads$gene %in% overlap_gene$gene)
head(reads_eoverlap_gene)
write.table(reads_eoverlap_gene, "Normallized_reads_SelpKO_AS_genes.csv", sep=",",row.names = F)

###############################pheatmap
##############################input normalized data with selp
data<-read.csv("Normallized_reads_SelpKO_AS_genes.csv",header = T)
head(data)
##############################input normalized data without selp
data<-read.csv("Normallized_reads_SelpKO_low_high_output_genes_no_selp.csv",header = T)
head(data)
#######
data<-data[! duplicated(data$gene),]
data<-subset(data, data$gene !="")
data.matrix(data, rownames.force = NA)
rownames(data)<-data$gene
head(data)

#reorder using clustering callback fuction  
callback = function(hc, data){
  sv = svd(t(data))$v[,1]
  dend = reorder(as.dendrogram(hc), wts = sv)
  as.hclust(dend)
}
#reorder using dendsort package 
callback = function(hc, ...){dendsort(hc)}



#colum annotation
annotation_col = data.frame(
  Samples = factor(rep(c("WT", "KO"), c(3,2))))
annotation_col
row.names(annotation_col) <- colnames(data[,-1])
annotation_col
#row annotation
annotation_row = data.frame(
  Genes = factor(rep(c("Maintenance", "Aging_up","Aging_down"), c(4,8,3))))
annotation_row
row.names(annotation_row) <- data$gene#set orders by the orderin "genes" list
annotation_row

#heatmap
pheatmap(data[,-1], scale="row",
         annotation_col = annotation_col,
         annotation_row = annotation_row,
         cluster_rows=T, cluster_cols=T,
         show_rownames = T,           #hide gene names
         #clustering_method_rows = "average",
         border=FALSE,
         cellwidth = 12, cellheight = 10,
         fontsize_row = 10,fontsize = 10,
         #legend_breaks = c(-1, 0, 1)
         #legend=F,
         #cutree_cols = 2,
         #clustering_callback = callback
         )





