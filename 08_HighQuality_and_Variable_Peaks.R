library(dplyr)
library(Seurat)
peak1<-readRDS("peak_cluster_proportions_20cluster.rds")
all_testmy<-function(x){
  x<-as.numeric(x)
  all(x<0.05)
}
a<-apply(peak1,1,all_testmy)
aa<-which(a==F)
peak1<-peak1[aa,]
write.table(rownames(peak1),"HIghQualityPeak_pro0.05.bed",col.names = F,row.names = F,sep = "\t",quote = F)

peak1<-as.matrix(peak1)
dge <- CreateSeuratObject(counts = peak1)
dge <- FindVariableFeatures(dge, selection.method = "vst", nfeatures = 50000)
VariablePeak<-VariableFeatures(dge)
head(VariablePeak)
write.table(VariablePeak,"VariablePeak_top50000.bed",col.names = F,row.names = F,sep = "\t",quote = F)
