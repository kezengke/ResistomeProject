#sorted boxplot for in/out patient Wilcoxon
rm(list = ls())

metaData<-read.csv("MetaBinnedWithInOutDonorTreatmentType.csv", header = T, row.names = 1)

pres<-rownames(metaData)[metaData$bins == "PRE"]
posts<-rownames(metaData)[metaData$bins == "D30"]

plotFun <- function(pres, posts, preR, postR, countsT, meta, tableType) {
  #sort wilcox-test results
  sortedPreR<-preR[order(preR$PreWilcox_p), , drop = F]
  sortedPostR<-postR[order(postR$PostWilcox_p), , drop = F]
  
  countsT<-countsT[, intersect(colnames(countsT), rownames(meta)), drop = F]
  meta<-meta[intersect(colnames(countsT), rownames(meta)), , drop = F]
  
  patients<-unique(meta$ID)
  
  metaPost<-meta[posts, , drop = F]
  
  metaPre<-meta[pres, , drop = F]
  metaPre<-metaPre[metaPre$ID %in% metaPost$ID, , drop = F]
  
  PreT<-countsT[, rownames(metaPre), drop = F]
  PostT<-countsT[, rownames(metaPost), drop = F]
  
  #sort counts tables based on ttest results
  PreT<-PreT[rownames(sortedPreR), rownames(metaPre), drop = F]
  PostT<-PostT[rownames(sortedPostR), rownames(metaPost), drop = F]
  
  #pre plots
  pdf(paste0("Plots/", tableType, "PreInOutPatientSortedBoxPlots(Wilcoxon).pdf"), width=14, height=7)
  par(mfrow=c(2, 4))
  par(mar=c(5,6,4,1)+.1)  
  
  label = c("In", "Out")
  for (i in 1:nrow(countsT)){
    boxplot(unlist(PreT[i,])~metaPre$ptInOut, outline = F,
            ylab=rownames(PreT)[i], xlab="Pre Patient", 
            names = label,
            cex.lab = 1.5, cex.axis = 1, cex.main = 1.8, 
            col=c( "darkorange", "cornflowerblue"))
    stripchart(unlist(PreT[i,])~metaPre$ptInOut, method="jitter", 
               vertical=T, pch=19, add=T)
    
    mtext(paste("Wilcox P-value:",signif(preR[rownames(PreT[i,]),1],digits = 3)),3, 0.9, cex = 0.6, adj = 0)
    mtext(paste("Adj. Wilcox P-value:",signif(preR[rownames(PreT[i,]),2],digits = 3)),3, 0.9, cex = 0.6, adj = 1)
  }
  dev.off()
  
  #post plots
  pdf(paste0("Plots/", tableType, "PostInOutPatientSortedBoxPlots(Wilcoxon).pdf"), width=14, height=7)
  par(mfrow=c(2, 4))
  par(mar=c(5,6,4,1)+.1)  
  
  label = c("In", "Out")
  for (i in 1:nrow(countsT)){
    boxplot(unlist(PostT[i,])~metaPost$ptInOut, outline = F,
            ylab=rownames(PostT)[i], xlab="Post Patient", 
            names = label,
            cex.lab = 1.5, cex.axis = 1, cex.main = 1.8, 
            col=c( "darkorange", "cornflowerblue"))
    stripchart(unlist(PostT[i,])~metaPost$ptInOut, method="jitter", 
               vertical=T, pch=19, add=T)
    
    mtext(paste("Wilcox P-value:",signif(preR[rownames(PostT[i,]),1],digits = 3)),3, 0.9, cex = 0.6, adj = 0)
    mtext(paste("Adj. Wilcox P-value:",signif(postR[rownames(PostT[i,]),2],digits = 3)),3, 0.9, cex = 0.6, adj = 1)
  }
  dev.off()
}

#wilcox results
speciesPreResults<-read.csv("Wilcox/brackenPreInOut.csv", row.names = 1)
speciesPostResults<-read.csv("Wilcox/brackenPostInOut.csv", row.names = 1)
genusPreResults<-read.csv("Wilcox/genusPreInOut.csv", row.names = 1)
genusPostResults<-read.csv("Wilcox/genusPostInOut.csv", row.names = 1)
phylumPreResults<-read.csv("Wilcox/phylumPreInOut.csv", row.names = 1)
phylumPostResults<-read.csv("Wilcox/phylumPostInOut.csv", row.names = 1)

amrPreResults<-read.csv("Wilcox/amrPreInOut.csv", row.names = 1)
amrPostResults<-read.csv("Wilcox/amrPostInOut.csv", row.names = 1)
rgiPreResults<-read.csv("Wilcox/rgiPreInOut.csv", row.names = 1)
rgiPostResults<-read.csv("Wilcox/rgiPostInOut.csv", row.names = 1)
vsearchPreResults<-read.csv("Wilcox/vsearchPreInOut.csv", row.names = 1)
vsearchPostResults<-read.csv("Wilcox/vsearchPostInOut.csv", row.names = 1)

pathwayPreResults<-read.csv("Wilcox/pathwayPreInOut.csv", row.names = 1)
pathwayPostResults<-read.csv("Wilcox/pathwayPostInOut.csv", row.names = 1)

#counts tables
speciesT<-read.csv("CountsTables/brackenFiltered.csv", header = T, row.names = 1, check.names = F)
genusT<-read.csv("CountsTables/genusFiltered.csv", header = T, row.names = 1, check.names = F)
phylumT<-read.csv("CountsTables/phylumFiltered.csv", header = T, row.names = 1, check.names = F)
amrT<-read.csv("CountsTables/amrFiltered.csv", header = T, row.names = 1, check.names = F)
rgiT<-read.csv("CountsTables/rgiFiltered.csv", header = T, row.names = 1, check.names = F)
vsearchT<-read.csv("CountsTables/vsearchFiltered.csv", header = T, row.names = 1, check.names = F)
pathwayT<-read.csv("CountsTables/pathFiltered.csv", header = T, row.names = 1, check.names = F)

plotFun(pres, posts, speciesPreResults, speciesPostResults, speciesT, metaData, "Species")
plotFun(pres, posts, genusPreResults, genusPostResults, genusT, metaData, "Genus")
plotFun(pres, posts, phylumPreResults, phylumPostResults, phylumT, metaData, "Phylum")

plotFun(pres, posts, amrPreResults, amrPostResults, amrT, metaData, "AMR")
plotFun(pres, posts, rgiPreResults, rgiPostResults, rgiT, metaData, "RGI")
plotFun(pres, posts, vsearchPreResults, vsearchPostResults, vsearchT, metaData, "Vsearch")

plotFun(pres, posts, pathwayPreResults, pathwayPostResults, pathwayT, metaData, "Pathway")
