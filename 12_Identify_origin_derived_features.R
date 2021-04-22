##identify feature genes for pRCC
library("edgeR")
count<- read.csv("pRCC_cicero_gene_activities.csv")
info<-read.csv("Group information of 34 pRCC samples.csv",header = T)

group <- factor(info$cluster)
y <- DGEList(counts=count, group=group)
y <- calcNormFactors(y,method="TMM")
design<-NULL
estimateGLMCommonDisp(y, design,verbose = TRUE)
fit1 <- glmFit(y,design,dispersion=0.1291805)
lrt1 <- glmLRT(fit1)

lrt_m<-as.matrix(lrt1$table)
diff_matrix<-cbind(count,lrt_m)
diff_matrix<-as.data.frame(diff_matrix)
down<-subset(diff_matrix,diff_matrix$logFC < (-1) & diff_matrix$PValue<0.01 & diff_matrix$logCPM>5)
up<-subset(diff_matrix,diff_matrix$logFC > 1 & diff_matrix$PValue<0.01 & diff_matrix$logCPM>5)

down<-down[,35:38]
up<-up[,35:38]
write.csv(diff_matrix,"diff_matrix_pRCC.csv",row.names = T)
write.csv(down,"pRCCa2_marker.csv")
write.csv(up,"pRCCa1_marker.csv")

##identify feature genes for normal kidney cell types
activity<-readRDS("CD_ICA_PT_pro_diff.rds")
info<-read.csv("cluster_info_new_11group.csv",header = T)
head(info)
library(dplyr)
library(Seurat)
dge <- CreateSeuratObject(counts = activity, min.cells = 10)
dge <- NormalizeData(dge, normalization.method = "LogNormalize", scale.factor = 10000)
dge <- FindVariableFeatures(dge, selection.method = "vst", nfeatures = 5000)
all.genes <- rownames(dge)
dge <- ScaleData(dge, features = all.genes)
index<-c(5,6,9)
a<-which(info$cluster %in% index)
info1<-info[a,]
a<-which(info1$cluster=="9")
info1$cluster[a]<-"5"
cluster<-as.character(info1$cluster)
names(cluster)<-info1$barcode
dge[["celltype"]] <- cluster
cluster<-as.data.frame(cluster)
dge<-AddMetaData(dge,metadata = cluster)
cluster<-dge[["celltype"]]
cluster1<-as.character(cluster$celltype)
names(cluster1)<-rownames(cluster)
cluster1<-as.factor(cluster1)
dge@active.ident<-cluster1
dge.markers <- FindAllMarkers(dge, only.pos = TRUE)
allmakers<-dge.markers[dge.markers$pct.1>0.5 & dge.markers$avg_logFC>0.5,]
DT_marker<-allmakers[allmakers$cluster=="6",]
DT_marker$fc<-DT_marker$pct.1/DT_marker$pct.2
write.csv(DT_marker,"CD_PC_marker.csv")
PT_marker<-allmakers[allmakers$cluster=="5",]
PT_marker$fc<-PT_marker$pct.1/PT_marker$pct.2
write.csv(PT_marker,"PT_marker.csv")
saveRDS(dge,"CD_PC_PT_pro_geneactivityDiffAnalysis.rds")
write.csv(dge.markers,"CD_PC_PT_pro_geneactivity_AlldgeMarkers.csv")

##identify origin-derived features
CDPC_marker<-read.csv("CD_PC_marker.csv",header = T)
PT_marker<-read.csv("PT_marker.csv",header = T)
pRCCa1_marker<-read.csv("pRCCa1_marker.csv",header = T)
pRCCa2_marker<-read.csv("pRCCa2_marker1.csv",header = T)
pRCCa1_marker_in_PT<-pRCCa1_marker[pRCCa1_marker$X %in% PT_marker$X,]
pRCCa2_marker_in_PC<-pRCCa2_marker[pRCCa2_marker$X %in% CDPC_marker$X,]
write.csv(pRCCa1_marker_in_PT,"pRCCa1_marker_in_PT.csv")
write.csv(pRCCa2_marker_in_PC,"pRCCa2_marker_in_PC.csv")
