library(SnapATAC);
###Step 1. Create snap object
file.list = c("ZY.snap", "ZY-2.snap");
sample.list = c("ZY-1", "ZY-2");
x.sp.ls = lapply(seq(file.list), function(i){
  x.sp = createSnap(file=file.list[i], sample=sample.list[i]);
  x.sp
})
names(x.sp.ls) = sample.list;
sample.list
###Step 2. Select barcode
barcode.file.list = c("ZY-1_barcode.txt", "ZY-2_barcode.txt");
barcode.list = lapply(barcode.file.list, function(file){
  read.table(file)[,1];
})
x.sp.list = lapply(seq(x.sp.ls), function(i){
  x.sp = x.sp.ls[[i]];
  x.sp  = x.sp[x.sp@barcode %in% barcode.list[[i]],];
})
names(x.sp.list) = sample.list;
###Step 3. Add cell-by-bin matrix
x.sp.list = lapply(seq(x.sp.list), function(i){
  x.sp = addBmatToSnap(x.sp.list[[i]], bin.size=5000);
  x.sp
})
x.sp.list
dim(x.sp.list[[2]]@bmat)
head(x.sp.list[[2]]@bmat)[,1:5]
###Step 4. Combine snap objects
bin.shared = Reduce(intersect, lapply(x.sp.list, function(x.sp) x.sp@feature$name));
x.sp.list <- lapply(x.sp.list, function(x.sp){
  idy = match(bin.shared, x.sp@feature$name);
  x.sp[,idy, mat="bmat"];
})
dim(x.sp@bmat)
dim(as.data.frame(x.sp@feature))
x.sp@feature[50000:50005]
head(x.sp@feature)
x.sp = Reduce(snapRbind, x.sp.list);
rm(x.sp.list); # free memory
gc();
table(x.sp@sample);
###Step 5. Binarize matrix
x.sp = makeBinary(x.sp, mat="bmat");
###Step 6. Filter bins
library(GenomicRanges);
black_list = read.table("hg38.blacklist.bed.gz");
head(black_list)
black_list.gr = GRanges(
  black_list[,1], 
  IRanges(black_list[,2], black_list[,3])
);
idy1 = queryHits(findOverlaps(x.sp@feature, black_list.gr));
idy2 = grep("chrM|random|chrX|chrY|chrUn|GL|alt", x.sp@feature);
idy = unique(c(idy1, idy2));
if(length(idy) > 0){x.sp = x.sp[,-idy, mat="bmat"]};
x.sp
#remove the top 5% bins
bin.cov = log10(Matrix::colSums(x.sp@bmat)+1);
dim(x.sp@bmat)
head(x.sp@bmat)[,1:5]
range(bin.cov)
head(bin.cov)
bin.cutoff = quantile(bin.cov[bin.cov > 0], 0.95);
idy = which(bin.cov <= bin.cutoff & bin.cov > 0);
x.sp = x.sp[, idy, mat="bmat"];
x.sp
###Step 7. Dimensionality Reduction
#LSA-logTF modified by Andrew Hill
system.time({
  x.lsa.logTF.sp = runLSA(
    obj=x.sp, 
    input.mat="bmat", 
    pc.num=50, 
    logTF=TRUE,
    scale.factor=1e+05, 
    min.cell=10, 
    seed.use=10
  );
})


###Step 8. Remove batch effect
library(harmony);
x.lsa.logTF.after.sp = runHarmony(
  obj=x.lsa.logTF.sp, 
  pca.dims=2:50, 
  meta_data=x.sp@sample,
  weight.by.sd=T
)

###Step 9. KNN Graph Construction (SnapATAC)
x.lsa.logTF.after.sp = runKNN(
  obj=x.lsa.logTF.after.sp,
  pca.dims=2:40,
  weight.by.sd=TRUE,
  k=15
)


###Step 10. Clustering (SnapATAC)
x.lsa.logTF.after.sp = runCluster(
  obj=x.lsa.logTF.after.sp,
  tmp.folder=tempdir(),
  louvain.lib="R-igraph",
  seed.use=10
)
x.lsa.logTF.after.sp@metaData$cluster = x.lsa.logTF.after.sp@cluster

###Step 11. Visualization
#lsa.logTF
x.lsa.logTF.sp = runViz(
  obj=x.lsa.logTF.sp, 
  tmp.folder=tempdir(),
  dims=2,
  pca.dims=2:40, 
  method="Rtsne",
  seed.use=10
);
x.lsa.logTF.after.sp = runViz(
  obj=x.lsa.logTF.after.sp, 
  tmp.folder=tempdir(),
  dims=2,
  pca.dims=2:40, 
  method="Rtsne",
  seed.use=10
);
pdf("ZY_lsa.logTF.pdf")
par(mfrow = c(2, 3));
plotViz(
  obj=x.lsa.logTF.sp,
  method="tsne", 
  main="Before Harmony",
  point.color="sample", 
  point.size=0.1, 
  text.add= FALSE,
  down.sample=10000,
  legend.add=TRUE
);
plotViz(
  obj=x.lsa.logTF.after.sp,
  method="tsne", 
  main="After Harmony",
  point.color="sample", 
  point.size=0.1, 
  text.add=FALSE,
  down.sample=10000,
  legend.add=TRUE
);
plotViz(
  obj=x.lsa.logTF.after.sp,
  method="tsne", 
  main="Cluster",
  point.color="cluster", 
  point.size=0.1, 
  text.add=TRUE,
  text.size=1,
  text.color="black",
  text.halo.add=TRUE,
  text.halo.color="white",
  text.halo.width=0.2,
  down.sample=10000,
  legend.add=FALSE
)
dev.off()

###cluster_info
b<-as.data.frame(x.lsa.logTF.after.sp@barcode)
c<-as.data.frame(x.lsa.logTF.after.sp@cluster)
cluster_info<-cbind(b,c)
colnames(cluster_info)<-c("barcode","cluster")
write.csv(cluster_info,"cluster_info_lsa.logTF.csv",row.names = F)

saveRDS(x.lsa.logTF.after.sp,"x.lsa.logTF.after.sp.rds")
saveRDS(x.lsa.logTF.sp,"x.lsa.logTF.sp.rds")
saveRDS(x.sp,"x.sp.rds")
