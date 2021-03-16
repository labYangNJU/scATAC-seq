## this script was downloaded from https://github.com/shendurelab/fly-atac with a little change
library(argparse)
library(dplyr)
library(methods)

insert_sizes = "ZY_Insertsize.txt" ##ZY-2_insertsize.txt


# FFT-based nucleosome pattern scoring
get_banding_score = function(cell_subset) {
  size_range=0:1000
  
  cell_subset.cleaned = cell_subset[, c("cell", "insert_size", "read_count")]
  missing_rows = do.call(rbind, lapply(size_range[! size_range %in% cell_subset.cleaned$insert_size], function(x) { data.frame(cell=cell_subset.cleaned[1, "cell"], insert_size=x, read_count=0)}))
  
  
  cell_subset.cleaned = rbind(cell_subset.cleaned, missing_rows) %>% arrange(insert_size)
  
  periodogram = spec.pgram(cell_subset.cleaned$read_count / max(cell_subset.cleaned$read_count), pad=0.3, tap=0.5, span=2, plot=F, fast=T)
  
  periodogram$freq = 1 / periodogram$freq
  
  #banding_score = sum(periodogram$spec[periodogram$freq >= 100 & periodogram$freq <= 300])
  banding_score = sum(periodogram$spec[periodogram$freq >= 100 & periodogram$freq <= 300])/sum(periodogram$spec)
  
  return(data.frame("cell"=cell_subset[1, "cell"], "banding_score"=banding_score))
}

# Read data
insert_sizes = read.delim(insert_sizes)
head(insert_sizes)

insert_sizes$cell = factor(insert_sizes$cell)

## Split each cell into it's own dataset and then calculate a score for each
insert_sizes = insert_sizes %>% arrange(cell)
indices = 1:nrow(insert_sizes)
indices_by_cell = split(indices, insert_sizes$cell)
nucleosome_banding_scores = do.call(rbind, lapply(1:length(indices_by_cell), function(i) {
  print(i);
  index_set = indices_by_cell[[i]];
  data = insert_sizes[index_set, ]
  get_banding_score(data)
}))
head(nucleosome_banding_scores)
nucleosome_banding_scores<-nucleosome_banding_scores[nucleosome_banding_scores$cell %in% select_count_tss$barcode,]
colnames(nucleosome_banding_scores)<-c("barcode","banding_score")
select_count_tss<-read.csv("ZY_info_readcount3.4_tssRatio0.15.csv",header = T)
select_count_tss_insertsize<-merge(select_count_tss,nucleosome_banding_scores,by="barcode")
dim(select_count_tss_insertsize)

range(log10(select_count_tss_insertsize$banding_score))
library(mclust)
#### Determine insert size cutoff
cellcall3 = Mclust(data.frame(log10(select_count_tss_insertsize$banding_score)),G=3)
summary(cellcall3)
banding_scores_cutoff = min(log10(select_count_tss_insertsize[which(cellcall3$classification == 2),5]))






