ex<-read.delim("TCGA-KIRP.htseq_fpkm.tsv",header = T,row.names = 1)
ex<-as.data.frame(ex)
ex$id<-rownames(ex)
ex_nopoint <- ex %>% 
  tidyr::separate(id,into = c("id"),sep="\\.")


##
group<-read.delim("KIRP_255samples_markerGene_FPKM_scale_for_cluster_K_G2.kgg",header = T)
CIMP<-read.csv("KIRP_CIMP_clinical_info.csv",header = T)
CIMP$ID<-paste0(CIMP$ID,"A")
a<-match(CIMP$ID,group$sample)
group_CIMP<-cbind(group[a,],CIMP)

##create matrix
gene_anno <- rtracklayer::readGFF("Homo_sapiens.GRCh38.100.gtf")
gene_anno<-gene_anno[gene_anno$type=="gene",]
gene_anno<-gene_anno[gene_anno$gene_biotype=="protein_coding",]
a<-which(ex_nopoint$id %in% gene_anno$gene_id)
ex_nopoint<-ex_nopoint[a,]
rownames(ex_nopoint)<-ex_nopoint$id
ex_nopoint1<-ex_nopoint[,-322]

group1_sample<-as.character(group$sample[1:41])
a<-which(group1_sample %in% CIMP$ID)
NotCIMP_sample<-group1_sample[-a]
CIMP_sample<-group1_sample[a]
NotCIMP_ex<-ex_nopoint1[,colnames(ex_nopoint1) %in% NotCIMP_sample]
CIMP_ex<-ex_nopoint1[,colnames(ex_nopoint1) %in% CIMP_sample]
pRCCa1_sample<-as.character(group$sample[42:255])
pRCCa1_Like_ex<-ex_nopoint1[,colnames(ex_nopoint1) %in% pRCCa1_sample]

normal<-read.csv("normal_info.csv",header = T)
normal$id<-gsub("-",".",normal$id)
normal_ex<-ex_nopoint1[,colnames(ex_nopoint1) %in% normal$id]
fpkm_ex_full<-cbind(CIMP_ex,NotCIMP_ex,pRCCa1_Like_ex,normal_ex)


##trajectory
pRCCa1_Like_ex1<-as.matrix(2^pRCCa1_Like_ex-1)
pRCCa1_Like_deviation<-apply(pRCCa1_Like_ex1,1,sd)
pRCCa1_Like_mean<-apply(pRCCa1_Like_ex1,1,mean)
a<-which(pRCCa1_Like_deviation<10 & pRCCa1_Like_mean>1)
pRCCa1_Like_ex1<-pRCCa1_Like_ex1[a,]

CIMP_ex1<-as.matrix(2^CIMP_ex-1)
CIMP_deviation<-apply(CIMP_ex1,1,sd)
CIMP_mean<-apply(CIMP_ex1,1,mean)
a<-which(CIMP_deviation<10 & CIMP_mean>1)
CIMP_ex1<-CIMP_ex1[a,]

NotCIMP_ex1<-as.matrix(2^NotCIMP_ex-1)
NotCIMP_deviation<-apply(NotCIMP_ex1,1,sd)
NotCIMP_mean<-apply(NotCIMP_ex1,1,mean)
a<-which(NotCIMP_deviation<10 & NotCIMP_mean>1)
NotCIMP_ex1<-NotCIMP_ex1[a,]

normal_ex1<-as.matrix(2^normal_ex-1)
normal_deviation<-apply(normal_ex1,1,sd)
normal_mean<-apply(normal_ex1,1,mean)
a<-which(normal_deviation<10 & normal_mean>1)
normal_ex1<-normal_ex1[a,]

use_gene<-unique(c(rownames(CIMP_ex1),rownames(NotCIMP_ex1),rownames(pRCCa1_Like_ex1),rownames(normal_ex1)))
trajectory_ex<-count[rownames(count) %in% use_gene,]

##
sample_sheet<-as.data.frame(colnames(trajectory_ex))
sample_sheet$celltype<-c(rep("CIMP",8),rep("NotCIMP",33),rep("pRCCa1_Like",214),rep("Normal",32))
sample_sheet1<-as.data.frame(sample_sheet[,-1])
rownames(sample_sheet1)<-colnames(trajectory_ex)
colnames(sample_sheet1)<-"celltype"
a<-match(rownames(trajectory_ex),gene_anno$gene_id)
gene_annotation<-as.data.frame(gene_anno$gene_name[a])
rownames(gene_annotation)<-rownames(trajectory_ex)
colnames(gene_annotation)<-"gene_short_name"

library(monocle)
pd <- new("AnnotatedDataFrame", data = sample_sheet1)
fd <- new("AnnotatedDataFrame", data = gene_annotation)
cds <- newCellDataSet(as.matrix(trajectory_ex),phenoData = pd, featureData = fd,expressionFamily=negbinomial.size())
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
cds <- detectGenes(cds, min_expr = 0.1)
expressed_genes <- row.names(subset(fData(cds),num_cells_expressed >= 4))
##celltype
pData(cds)$CellType<-sample_sheet1$celltype
pData(cds)$CellType2<-"pRCCa1_Like"
pData(cds)$CellType2[1:41]<-"pRCCa2_Like"
pData(cds)$CellType2[256:287]<-"Normal"
#pie plot
pie <- ggplot(pData(cds),aes(x = factor(1), fill = factor(CellType))) + geom_bar(width = 1)
pie + coord_polar(theta = "y") +theme(axis.title.x = element_blank(), axis.title.y = element_blank())

#cluster cells
cds <- reduceDimension(cds, max_components = 2, num_dim = 4,
                       norm_method = 'log',
                       reduction_method = 'tSNE',
                       residualModelFormulaStr = "~num_genes_expressed",
                       verbose = T)
cds <- clusterCells(cds, num_clusters = 4)
plot_cell_clusters(cds, 1, 2,3,4, color = "CellType")

##differential genes
diff_test_res <- differentialGeneTest(cds[expressed_genes,],fullModelFormulaStr = "~CellType")

#ordering_genes <- row.names (diff_test_res)
ordering_genes <- row.names (diff_test_res[diff_test_res$qval < 0.00003,])

#marks genes that will be used for clustering in subsequent calls
cds <- setOrderingFilter(cds, ordering_genes)
#reduce data dimensionality
cds <- reduceDimension(cds, max_components = 2,method = 'DDRTree')
#order cells along the trajectory 
cds <- orderCells(cds)
plot_cell_trajectory(cds, color_by = "CellType",show_branch_points = F)
plot_cell_trajectory(cds, color_by = "State")
