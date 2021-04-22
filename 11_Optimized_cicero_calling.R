gene_anno <- rtracklayer::readGFF("Homo_sapiens.GRCh38.100.gtf")
gene_anno$chromosome <- paste0("chr", gene_anno$seqid)
gene_anno$gene <- gene_anno$gene_id
gene_anno$transcript <- gene_anno$transcript_id
gene_anno$symbol <- gene_anno$gene_name
gene_anno<-gene_anno[,c(27,3:5,7,11,13)]
gene_anno<-gene_anno[gene_anno$type=="gene",]
gene_anno<-gene_anno[gene_anno$gene_biotype=="protein_coding",]
gene_anno_name<-gene_anno
pos<-which(gene_anno_name$strand=="+")
gene_anno_name$start_new<-gene_anno_name$start
gene_anno_name$start_new[pos]<-gene_anno_name$start_new[pos]-200
gene_anno_name$end_new<-gene_anno_name$start
gene_anno_name$end_new[pos]<-gene_anno_name$end_new[pos]+100
neg<-which(gene_anno_name$strand=="-")
gene_anno_name$start_new[neg]<-gene_anno_name$end[neg]-100
gene_anno_name$end_new[neg]<-gene_anno_name$end[neg]+200
gene_anno_name$peakname<-paste(gene_anno_name$chromosome,gene_anno_name$start_new,gene_anno_name$end_new,sep = "_")
gene<-gene_anno_name[,c(1,3,4)]
colnames(gene)[1:3]<-c("V1","V2","V3")
library(GenomicRanges)
gene_bed<-with(gene,GRanges(V1,IRanges(V2,V3)))
peakinfo <- read.table("HIghQualityPeak_pro0.05.bed",header = F)
peakinfo<-peakinfo[,1:3]
names(peakinfo) <- c("chr", "bp1", "bp2")
peakinfo$site_name <- paste(peakinfo$chr, peakinfo$bp1, peakinfo$bp2, sep="_")
rownames(peakinfo)<-peakinfo$site_name
ZYpeak_bed<-with(peakinfo,GRanges(chr,IRanges(bp1,bp2)))
o<-findOverlaps(gene_bed,ZYpeak_bed)
a<-unique(as.numeric(queryHits(o)))
gene_overlap<-gene_anno_name[a,]
gene_overlap_peak<-gene_overlap[,c(1,3,4)]
colnames(gene_overlap_peak)[1:3]<-c("V1","V2","V3")
gene_overlap_peak_bed<-with(gene_overlap_peak,GRanges(V1,IRanges(V2,V3)))
##construct genebody peak info
geneinfo<-gene_overlap[,c(1,8:10)]
a<-duplicated(geneinfo$peakname)
b<-which(a==T)
geneinfo$end_new[b]<-geneinfo$end_new[b]+1
geneinfo$peakname<-paste(geneinfo$chromosome, geneinfo$start_new, geneinfo$end_new, sep="_")
rownames(geneinfo)<-geneinfo$peakname
colnames(geneinfo)<-c("chr","bp1","bp2","site_name")
human.hg38.genome<-read.delim("hg38.chrom.sizes",header = F)
library(monocle)
library(cicero)
peakinfo<-rbind(peakinfo,geneinfo)

###cluster2 for example
indata <- Matrix::readMM("cluster2/matrix.mtx")
indata@x[indata@x > 0] <- 1
cellinfo <- read.table("cluster2/barcode.txt")
row.names(cellinfo) <- cellinfo$V1
names(cellinfo) <- "cells"
colnames(indata)<-cellinfo$cells

colnames(indata) <- row.names(cellinfo)
gene_matrix<-matrix(0,ncol=dim(cellinfo)[1],nrow=16440)
for(i in 1:16440){
  gene_tmp<-gene_overlap_peak_bed[i,]
  o<-findOverlaps(gene_tmp,ZYpeak_bed)
  hit<-as.numeric(subjectHits(o))
  if(length(hit)==1){
    hit_matrix<-t(as.matrix(indata[hit,]))
    gene_matrix[i,]<-hit_matrix}
  else{
    hit_matrix<-indata[hit,]
    dim(hit_matrix)
    hit_matrix<-as.matrix(hit_matrix)
    count<-t(as.matrix(colSums(hit_matrix)))
    length(count)
    gene_matrix[i,]<-count}
}

rownames(gene_matrix)<-geneinfo$site_name
colnames(gene_matrix)<-colnames(indata)

gene_matrix<-Matrix(gene_matrix,sparse = T)
indata<-rbind(indata,gene_matrix)
rownames(indata)<-peakinfo$site_name

fd <- methods::new("AnnotatedDataFrame", data = peakinfo)
pd <- methods::new("AnnotatedDataFrame", data = cellinfo)
input_cds <-  suppressWarnings(newCellDataSet(indata,
                                              phenoData = pd,
                                              featureData = fd,
                                              expressionFamily=VGAM::binomialff(),
                                              lowerDetectionLimit=0))
input_cds@expressionFamily@vfamily <- "binomialff"
input_cds <- monocle::detectGenes(input_cds)
input_cds <- input_cds[Matrix::rowSums(exprs(input_cds)) != 0,]
set.seed(2017)
input_cds <- detectGenes(input_cds)
input_cds <- estimateSizeFactors(input_cds)
input_cds <- reduceDimension(input_cds, max_components = 2, num_dim=6,
                             reduction_method = 'tSNE', norm_method = "none")
tsne_coords <- t(reducedDimA(input_cds))
row.names(tsne_coords) <- row.names(pData(input_cds))
cicero_cds <- make_cicero_cds(input_cds, reduced_coordinates = tsne_coords,k=30)
##Call Cicero
conns <- run_cicero(cicero_cds, human.hg38.genome)
conns1<-conns[conns$coaccess>0,]
index<-which(is.na(conns1))
head(index)
conns1=conns1[-index,]
saveRDS(cicero_cds,"cluster2/cluster2_cicero_cds.rds")
saveRDS(input_cds,"cluster2/cluster2_input_cds.rds")
write.csv(conns1,"cluster2/cluster2_conns.csv",row.names = F)
conns2<-conns1[conns1$coaccess>0.1,]
write.csv(conns2,"cluster2/cluster2_conns_cutoff0.1.csv",row.names = F)


##calculate gene activity matrix
input_cds<-readRDS("cluster2/cluster2_input_cds.rds")
conns<-read.csv("cluster2/cluster2_conns_cutoff0.1.csv",header = T)
a<-match(fData(input_cds)$site_name,geneinfo$site_name)
fData(input_cds)$gene<-gene_overlap$gene_name[a]
unnorm_ga <- build_gene_activity_matrix(input_cds, conns)
unnorm_ga <- unnorm_ga[!Matrix::rowSums(unnorm_ga) == 0,]
write.csv(as.matrix(unnorm_ga),"cluster2/cluster2_cicero_gene_activities.csv")
