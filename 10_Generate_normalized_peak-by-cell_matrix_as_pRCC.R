library(Rsamtools)
library(GenomicRanges)
p2 <- ScanBamParam(what=c("rname", "strand", "pos", "qwidth"))
cluster1<-scanBam("ZY_cluster1_logTF.bam",param=p2,maxMemory=5,package="Rsamtools")

##corrected for the Tn5 offset
strand<-cluster1[[1]]$strand
index<-which(as.character(strand)=="+")
pos<-cluster1[[1]]$pos
pos[index]<-pos[index]+4
cluster1[[1]]$start[index]<-pos[index]-cluster1[[1]]$qwidth[index]
cluster1[[1]]$end[index]<-pos[index]
index<-which(as.character(strand)=="-")
pos[index]<-pos[index]-5
cluster1[[1]]$start[index]<-pos[index]
cluster1[[1]]$end[index]<-pos[index]+cluster1[[1]]$qwidth[index]
cluster1[[1]]$pos<-pos
b<-which(cluster1[[1]]$rname=="Y")

fragement<-cbind(cluster1[[1]]$rname,cluster1[[1]]$start,cluster1[[1]]$end)
fragement<-as.data.frame(fragement)
fragement$V1<-paste0("chr",fragement$V1)
a<-which(fragement$V1=="chr23")
fragement$V1[a]="chrX"
a<-which(fragement$V1=="chr195")
fragement$V1[a]="chrY"
fragement_bed<-with(fragement, GRanges(V1, IRanges(V2, V3)))

##readin peak file
KIRP<-read.delim("../CancerPeak.bed",header = F)
KIRP_bed<-with(KIRP, GRanges(V1, IRanges(V2, V3)))

##countoverlap
tmp<-countOverlaps(KIRP_bed,fragement_bed)
count_matrix<-matrix(0,ncol = 20,nrow=dim(KIRP)[1])
count_matrix[,1]<-tmp

##repeat for the left clusters
name<-c("cluster1","cluster2","cluster3","cluster4","cluster5","cluster6","cluster7","cluster8","cluster9","cluster10","cluster11","cluster12","cluster13","cluster14","cluster15","cluster16","cluster17","cluster18","cluster19","cluster20")
colnames(count_matrix)<-name
rownames(count_matrix)<-KIRP$V4
dim(count_matrix)
head(count_matrix)[,1:5]
saveRDS(count_matrix,"ZY_count_matrix.rds")

##normalized by using edgeR’s cpm
library(edgeR)
count_matrix1<-cpm(count_matrix,log = TRUE, prior.count = 5)
##quantile normalization
biocLite("preprocessCore")
library(preprocessCore)
count_matrix2=normalize.quantiles(count_matrix1)
colnames(count_matrix2)<-colnames(count_matrix1)
write.table(count_matrix2,"ZY_normalizeMatrix.txt",row.names = F,sep = "\t",quote = F)

