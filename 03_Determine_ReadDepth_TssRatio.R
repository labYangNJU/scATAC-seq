info<-read.csv("ZY12_info_readcount_tssRatio_full.csv",header = T)
##method 1
library(mclust)
info10<-info[info$All_reads>10,]
hist(log10(info10$All_reads),breaks=60,col="mediumseagreen",main="Hist of read count per barcode",
     xlab="Number of Reads (log10)",las=1,xlim=c(0,6),ylab = "Num of barcode",cex.axis=0.6)
abline(v=log10(readCount_cutoff),col="red")
info10$All_reads<-as.numeric(info10$All_reads)
cellcall = Mclust(data.frame(log10(info10$All_reads)),G=2)
summary(cellcall)
a<-which(cellcall$classification == 2)
readCount_cutoff = min(info10[which(cellcall$classification == 2 & cellcall$uncertainty < 0.05),2])
#### Determine tss enrichment
hist(info10$tss_ratio,breaks=60,col="mediumseagreen",main="Hist of tss enrichment",
     xlab="tss enrichment",las=1,xlim=c(0,0.5),ylab = "Num of barcode",cex.axis=0.6)
abline(v=tss_enrichment_cut,col="red")
cellcall2 = Mclust(data.frame(info10$tss_ratio),G=2)
summary(cellcall2)
tss_enrichment_cut = min(info10[which(cellcall2$classification == 2 & cellcall2$uncertainty < 0.05),4])

##method 2 ----plot
info100<-info[info$All_reads>100,]
plot(log10(info100$All_reads),info100$tss_ratio,cex=0.1)
