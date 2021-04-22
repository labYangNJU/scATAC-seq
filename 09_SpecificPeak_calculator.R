## this script was downloaded from https://github.com/shendurelab/mouse-atac/tree/master/specificity with a little change
library(Matrix)
library(limma)
library(reshape2)
library(argparse)
library(methods)


makeprobsvec<-function(p){
  phat<-p/sum(p)
  phat[is.na(phat)] = 0
  phat
}

shannon.entropy <- function(p) {
  if (min(p) < 0 || sum(p) <=0)
    return(Inf)
  p.norm<-p[p>0]/sum(p)
  -sum( log2(p.norm)*p.norm)
}

JSdistVec <- function(p, q){
  JSdiv <- shannon.entropy((p + q)/2) - (shannon.entropy(p) + 
                                           shannon.entropy(q)) * 0.5
  JSdiv[is.infinite(JSdiv)] <- 1
  JSdiv[JSdiv < 0] <- 0
  JSdist <- sqrt(JSdiv)
  JSdist
}

specificity_scorer = function(normpropmat){
  marker_gene_specificities <- lapply(1:ncol(normpropmat), function(cell_type_i){
    perfect_specificity <- rep(0.0, ncol(normpropmat))
    perfect_specificity[cell_type_i] <- 1.0
    apply(normpropmat, 1, function(x) { 
      if (sum(x) > 0) 1 - JSdistVec(makeprobsvec(x), perfect_specificity)
      else 0
    })
  })
  return(do.call(cbind, marker_gene_specificities))
}

markerlistmaker = function(markering,daps){
  markerlist = list()
  betamarkerlist = list()
  nonmarkerlist = list()
  markerback = list()
  for(i in 1:ncol(markering)){
    commoners = intersect(rownames(daps[[match(colnames(markering)[i],names(daps))]]),rownames(markering))
    markerlist[[i]] = markering[match(commoners,rownames(markering)),i]
    betamarkerlist[[i]] = daps[[match(colnames(markering)[i],names(daps))]][match(commoners,rownames(daps[[match(colnames(markering)[i],names(daps))]])),1]
    nonmarkerlist[[i]] = markering[-match(commoners,rownames(markering)),i]
    names(markerlist)[i] = colnames(markering)[i]
    names(betamarkerlist)[i] = colnames(markering)[i]
    names(nonmarkerlist)[i] = colnames(markering)[i]
  }
  markerback$markerlist = markerlist
  markerback$betamarkerlist = betamarkerlist
  markerback$nonmarkerlist = nonmarkerlist
  return(markerback)
}

spec_floorer = function(trueo,nullo){
  nullpass = c()
  for(i in 1:length(trueo)){
    nullpass[i] = length(which(nullo >= trueo[i]))
  }
  nullfrac = nullpass/(nullpass + 1:length(nullpass))
  return(trueo[max(which(nullfrac <= 0.1))])
}

results_writer = function(markerlist,specfloor,daps,txtout,rdsout){
  wbmat = matrix(,,5)
  colnames(wbmat) = c("specificity_score","locusID","cluster","subset_cluster","cluster_name")
  for(i in 1:length(markerlist))
  {
    #print(i)
    specsites = markerlist[[i]][which(log10(markerlist[[i]]) > specfloor)]
    if(length(specsites) == 0){
      next
    }
    currmatout = matrix(names(specsites))
    currcluster = names(markerlist)[i]
    currmatout = cbind(currmatout,daps[[match(currcluster,names(daps))]][match(names(specsites),rownames(daps[[match(currcluster,names(daps))]])),])
    currmatout = cbind(currmatout,specsites)
    colnames(currmatout)[c(1,5)] = c("locusID","specificity_score")
    currmatout = currmatout[order(currmatout[,5],decreasing=T),]
    currmatout[,2:5] = signif(currmatout[,2:5],4)
    #print(dim(currmatout)[1])
    clusterer = strsplit2(currcluster,"[.]")
    clusternow = gsub("clusters_","",clusterer[1])
    subnow = gsub("cluster_","",clusterer[2])
    assigns = cbind(rep(clusternow,nrow(currmatout)),rep(subnow,nrow(currmatout)))
    assigns = cbind(assigns,rep(currcluster,nrow(currmatout)))
    currmatforwbmat = cbind(currmatout[,c(5,1)],assigns)
    rownames(currmatforwbmat) = NULL
    colnames(currmatforwbmat) = c("specificity_score","locusID","cluster","subset_cluster","cluster_name")
    wbmat = rbind(wbmat,currmatforwbmat)
  }
  wbmat = wbmat[-1,]
  write.table(wbmat,txtout,row.names=F,quote=F,sep="\t")
  saveRDS(wbmat,rdsout)
}

print("Loading files...")
propmat = readRDS("peak_cluster_proportions_11group_VariablePeak.rds")
#propmat<-1-propmat
depthmat = read.table("median_depth_logTF_20cluster_11group.txt")

print("Normalizing proportions...")
depthnorm = mean(depthmat[,2])/depthmat[,2]
logdepthnorm = mean(log10(depthmat[,2]))/log10(depthmat[,2])
propmatnormbylogdepth = t(t(propmat)*logdepthnorm)

print("Calculating specificity scores...")
marker_specificities_out = specificity_scorer(propmatnormbylogdepth)
dim(marker_specificities_out)
markerdup = marker_specificities_out^2
markering = markerdup * propmatnormbylogdepth
dim(markering)
rownames(markering) = rownames(propmat)
colnames(markering) = colnames(propmat)
markeringlong = melt(as.matrix(markering))
markeringlong = cbind(strsplit2(markeringlong[,2],"[.]"),markeringlong)
markeringlong = markeringlong[,c(5,3,1,2,4)]
colnames(markeringlong) = c("specificity_score","locusID","cluster","subset_cluster","cluster_name")

print("Writing all scores...")
markeringlong_s<-markeringlong[order(markeringlong$specificity_score,decreasing=TRUE),]
specific_top20000<-markeringlong_s[1:20000,]
write.csv(specific_top20000,"specific_top20000.csv",row.names = F)
