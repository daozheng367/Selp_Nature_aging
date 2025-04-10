
#install pakage edge

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("edgeR")

BiocManager::install("limma",force = TRUE)
BiocManager::install('EnhancedVolcano')

# ready to start

library(limma)
library(edgeR)
library(ggplot2)
library(ggrepel)
library(EnhancedVolcano)

mydata<-read.csv("EW210520_annt_donor.tsv", sep="\t", header=TRUE, row.names=1)
mydata
head(mydata)
row.names(mydata)
Gene<-mydata[,1]
anno<-data.frame(row.names(mydata), Gene, 
                 row.names=1)
anno
head(anno)

####end of annotation
#make DGElist, this is a special object in edgeR
df<-data.frame(mydata[,-1])
head(df)
y <- DGEList(df, remove.zeros = TRUE)     
head(y)

## Define samples
colnames(df)
express<-c(rep("O_L",4),rep("O_H",4))
#express<-c("O_L1","O_L2","O_L3","O_L4","O_H1","O_H2","O_H3","O_H4")
mouse<-c("M1","M2","M3","M4","M1","M2","M3","M4")
head(mouse)

express
#cells<-c(rep("Donor",8),rep("Transplant",6))
#condition<- as.factor(express) #choose either express or cells
condition<- as.factor(express)
condition
y$samples$condition <- relevel(condition, ref = "O_H") #any subgroup of selected

# Normalise (more precisely calculate normalizing factors)
y <- calcNormFactors(y,lib.size = y$samples$lib.size,method   = "TMM",p= 0.75)
head(y)
head(y$samples, length(express))
#make design for DGE
#design<-model.matrix(~0+ express)
design<-model.matrix(~0 + express+mouse)
#design<-model.matrix(~0+ condition +rep)
design

# Estimates the abundance-dispersion trend by Cox-Reid approximate profile likelihood.
y <- estimateGLMTrendedDisp(y,design)
plotBCV(y, main="GLMTrended")
#probably and
y <- estimateGLMTagwiseDisp(y,design,prior.df = 10)
plotBCV(y, main="GLMTagwise")
#OR
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
          #          contrast= makeContrasts(expressLowT- expressHighT, levels= design))
          #           contrast = makeContrasts(expressHigh-expressLow, levels = design))
                        contrast= makeContrasts(expressO_H-expressO_L, levels=design))
#OldVsYoungMinusRun <- glmLRT(fit, contrast = c(-1, 1))
# Identify which genes are significantly differentially expressed 
dtLowVsHigh <- decideTestsDGE(LowVsHigh, adjust.method = "BH")#, p.value = 0.05)

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
write.table(merged, "Volcano_Selp_GLMR_OL_VS_OH_input_paired.tsv", sep="\t", row.names=FALSE)



########Volcano plot

merged<-read.csv("Volcano_Selp_GLMR_OL_VS_OH_input_paired.tsv", sep="\t",header = TRUE)
head(merged)

keyvals <- ifelse(merged$logFC >0&merged$FDR <0.01, 'red',ifelse(merged$logFC < -0&merged$FDR <0.01, 'blue','black'))       
keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'red'] <- 'Up'
names(keyvals)[keyvals == 'black'] <- 'NS'
names(keyvals)[keyvals == 'blue'] <- 'Down'

EnhancedVolcano(merged,
                lab =merged$Gene ,
                x = 'logFC',
                y = 'PValue',
                #ylim = c(0,15),
                selectLab ="",
                title = "",
                subtitle = "",
                pCutoff = 0.01,
                pCutoffCol = 'FDR',#use FDR as cut off.Defult cut of is p value
                FCcutoff = 0,
                pointSize = 2.0,
                labSize = 3,
                colAlpha = 1,
                colCustom = keyvals,
                legendPosition = 'none',
                drawConnectors = TRUE,
                widthConnectors = 0.5,
                typeConnectors="open",
                )








