
########PSGL-1 treatment Volcano plot
library(ggplot2)
library(ggrepel)
library(EnhancedVolcano)
#load file
merged<-read.csv("DEGs_WTLG_WTVH_RNAseq.csv", sep=",", header=TRUE)
head(merged)

#volcano

keyvals <- ifelse(merged$log2FoldChange > 1&merged$padj <0.05, 'red',ifelse(merged$log2FoldChange < -1&merged$padj <0.05, 'blue','black'))


keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'red'] <- 'Up'
names(keyvals)[keyvals == 'black'] <- 'NS'
names(keyvals)[keyvals == 'blue'] <- 'Down'

EnhancedVolcano(merged,
                lab =merged$Gene,
                x = 'log2FoldChange',
                y = 'padj',
                #selectLab = c('Selp'),
                selectLab ="",#no gene lable
                #boxedLabels = TRUE,
                labFace = 'bold',
                labCol = 'black',
                #xlim = c(-30,30),
                #ylim = c(0,15),
                title = "",
                subtitle = "",
                pCutoff = 0.05,
                pCutoffCol = 'padj',#use FDR as cut off.Defult cut of is p value
                FCcutoff = 1,
                pointSize = 1,
                labSize = 5,
                axisLabSize = 15,
                colAlpha = 1,
                colCustom = keyvals,
                legendPosition = 'none',#+ coord_flip()
                #drawConnectors = TRUE,
                #widthConnectors = 0.5
                 )
                


