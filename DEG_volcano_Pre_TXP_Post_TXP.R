#Pre_TXP vs Post_TXP

library(limma)
library(edgeR)
library(ggplot2)
library(ggrepel)
library(EnhancedVolcano)
setwd("C:/PhD/Sequencing/Selp/DEG in R")
#Donor VS TXP
mydata<-read.csv("EW210520_annt.tsv", sep="\t", header=TRUE, row.names=1)
head(mydata)
row.names(mydata)
Gene<-mydata[,1]
anno<-data.frame(row.names(mydata), Gene, 
                 row.names=1)
head(anno)

####end of annotation
#make DGElist, this is a special object in edgeR
df<-data.frame(mydata[,-1])
head(df)
y <- DGEList(df, remove.zeros = TRUE)     
head(y)

## Define samples
colnames(df)
cells<-c(rep("Donor",8),rep("Transplant",6))

express<-c(rep("Low",4),rep("High",4),rep("Low",3),rep("High",3))
#express<-c("O_L1","O_L2","O_L3","O_L4","O_H1","O_H2","O_H3","O_H4")
mouse<-c("m1","m2","m3","m4","m1","m2","m3","m4","m4","m4","m4","m4","m4","m4")
condition<- as.factor(cells) #choose either express or cells
y$samples$condition <- relevel(condition, ref = "Transplant") #any subgroup of selected
par(mar=c(6,2,1,1))


# Normalise (more precisely calculate normalizing factors)
y <- calcNormFactors(y,lib.size = y$samples$lib.size,method   = "TMM",p= 0.75)
head(y)
head(y$samples, length(cells))
#make design for DGE
design<-model.matrix(~0 + cells+express+mouse)  #3design
#design<-model.matrix(~0+ condition +rep)
design


# Compute a robust estimate of the negative binomial dispersion parameter for each gene, with expression levels specified by a 
# log-linear model, using observation weights. 
y <- estimateGLMRobustDisp(y,design,prior.df = 10,verbose = TRUE)
plotBCV(y, main="GLMRobust")

# glmFit produces an object of class DGEGLM containing components 
# counts, samples, genes and abundance from y plus new components.
fit <- glmFit(y,design      = design,prior.count = 0.125)
fit$design
# Make DE lists again but not taken into account the run information
LowVsHigh <- glmLRT(fit, 
                        contrast= makeContrasts(cellsTransplant-cellsDonor, levels=design))
#OldVsYoungMinusRun <- glmLRT(fit, contrast = c(-1, 1))
# Identify which genes are significantly differentially expressed 
dtLowVsHigh <- decideTestsDGE(LowVsHigh, adjust.method = "BH", p.value = 0.05)

# This DE takes into account the run information (batch effect)
summary(dtLowVsHigh)
dtLowVsHigh
#write.table(dtLowVsHigh, "DE_Selp_GLMT_Donor_Transplant_Pvalue cutoff.tsv", sep="\t", row.names=FALSE)

# Extracts the top DE tags in a data frame for a given pair of groups, 
# ranked by p-value or absolute log-fold change.
DE <- topTags(LowVsHigh, adjust.method = "BH", n = Inf)$table
#add annotation

merged=merge(anno,DE,  by="row.names")
summary(merged)
head(merged)
write.table(merged, "Volcano_Selp_GLMR_donor VS TXP_input_3design.tsv", sep="\t", row.names=FALSE)

########Volcano plot
library(ggplot2)
library(ggrepel)
library(EnhancedVolcano)
merged<-read.csv("Volcano_Selp_GLMR_donor VS TXP_input_3design.tsv", sep="\t", header=TRUE)
head(merged)
keyvals <- ifelse(merged$logFC > 0&merged$FDR <0.001, 'red',ifelse(merged$logFC < 0&merged$FDR <0.001, 'blue','black'))
        
keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'red'] <- 'Up'
names(keyvals)[keyvals == 'black'] <- 'NS'
names(keyvals)[keyvals == 'blue'] <- 'Down'

EnhancedVolcano(merged,
                lab =merged$Gene ,
                x = 'logFC',
                y = 'PValue',
                selectLab = c('Selp'),
                boxedLabels = TRUE,
                labFace = 'bold',
                labCol = 'black',
                #ylim = c(0,15),
                title = "",
                subtitle = "",
                pCutoff = 0.001,
                pCutoffCol = 'FDR',#use FDR as cut off.Defult cut of is p value
                FCcutoff = 0,
                pointSize = 1,
                labSize = 3,
                colAlpha = 1,
                colCustom = keyvals,
                legendPosition = 'none')#+ coord_flip()
                #drawConnectors = TRUE,
                #widthConnectors = 0.5)
                

