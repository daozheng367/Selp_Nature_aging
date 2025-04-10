

library(limma)
library(edgeR)
library(ggplot2)
library(ggrepel)
library(EnhancedVolcano)

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
express<-c(rep("Selp_lo",4),rep("Selp_hi",4))#donor
express<-c(rep("Pre_TXP",8),rep("Post_TXP",6))#TXP

express

condition<- as.factor(express)
condition
y$samples$condition <- relevel(condition, ref = "Selp_hi") #any subgroup of selected
y$samples$condition <- relevel(condition, ref = "Post_TXP") #any subgroup of selected


#MDS

label<-c("lo_1","lo_2","lo_3","lo_4","hi_1","hi_2","hi_3","hi_4")#donor
label<-c(rep("Pre_TXP",8),rep("Post_TXP",6))#TXP

MDS<-plotMDS(y,col=c(rep("#7fbf7b",4),rep("#af8dc3",4)),
       # labels=1:14, #you can use any kind of labels here
        #labels=cells,
         labels=label,
       pch = 16,font=2,
       main = "")
MDSdata<-cbind(MDS$x, MDS$y)
colnames(MDSdata)[1:2]<-c("x","y")
df<-data.frame(MDSdata)

ggplot(df, aes(x,y,color=label)) +
  geom_point()+
  geom_text(aes(label=label),hjust = 0.7,nudge_y = 0.11,size=4,fontface = "bold")+#for points labels
  scale_shape_manual(values=c(16,17))+
  #guide_legend(title = Sample)+
  theme_bw()+
  theme(legend.position="none")#remove legend








