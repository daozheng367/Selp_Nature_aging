
# Go analysis was perform in DAVID website

install.packages("ENRICHPLOT")
BiocManager::install("DOSE",force = TRUE)
BiocManager::install("enrichplot",force = TRUE)
#BiocManager::install("clusterProfiler")

library(ggplot2)
library(dplyr)
library(stringi)


DAVID<-read.delim('TXP_down_GO_KEGG_PATHWAY.txt')

head(DAVID)

enrich_signif=DAVID[which(DAVID$FDR<0.05),]

enrich_signif=enrich_signif[,c(1:3,4,5)]

head(enrich_signif)

enrich_signif=data.frame(enrich_signif)

#head(enrich_signif)
##arrange order and sort

#rank top 10 by GeneRatio
enrich_signif_rank=arrange(enrich_signif,desc(enrich_signif$GeneRatio))[1:22,]#rank from biggest

#clean term names
enrich_signif_rank$Term<-stri_sub(enrich_signif_rank$Term,10,100)#start from number 12 character till 100, from 10 for KEGG, from 12 for GO terms
#head(enrich_signif_rank$Term)
#reoder based on GeneRatio
enrich_signif_rank$Term <- factor(enrich_signif_rank$Term, levels = enrich_signif_rank$Term[order(enrich_signif_rank$GeneRatio)])
#head(enrich_signif_rank$Term)
#ggplot
ggplot(enrich_signif_rank,aes(x=GeneRatio,y=Term))+
       labs(x = "Gene Ratio", y = "Down KEGG Terms", title = "") +
       geom_point(aes(color=PValue,size=Count))+
       scale_color_gradient(low='blue',high='red')+
       theme_bw()+
       theme(legend.position = c(1, 0),legend.justification = c(1,0))+#legend position
       theme(legend.key.size = unit(8, "pt"))+#legend size
       #labs(size=5)+
       theme(legend.background = element_blank())+#remove legend background
       theme(panel.grid.minor = element_blank(),
             panel.grid.major = element_blank(),
             axis.title.x = element_text(size = 12, face = "bold"),
             axis.text.x = element_text(size = 12, face = "bold"),
             axis.text.y = element_blank()#hide y lables
            )

