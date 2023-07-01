#sorted boxplot for in/out patient Ttest
rm(list = ls())

metaData<-read.csv("metaWithBins.csv", header = T, row.names = 1)

#wilcox results
brackenPreResults<-read.csv("Ttest/brackenPreInOut.csv", row.names = 1)
brackenPostResults<-read.csv("Ttest/brackenPostInOut.csv", row.names = 1)
amrPreResults<-read.csv("Ttest/amrPreInOut.csv", row.names = 1)
amrPostResults<-read.csv("Ttest/amrPostInOut.csv", row.names = 1)
rgiPreResults<-read.csv("Ttest/rgiPreInOut.csv", row.names = 1)
rgiPostResults<-read.csv("Ttest/rgiPostInOut.csv", row.names = 1)
vsearchPreResults<-read.csv("Ttest/vsearchPreInOut.csv", row.names = 1)
vsearchPostResults<-read.csv("Ttest/vsearchPostInOut.csv", row.names = 1)
genusPreResults<-read.csv("Ttest/genusPreInOut.csv", row.names = 1)
genusPostResults<-read.csv("Ttest/genusPostInOut.csv", row.names = 1)

#sort ttest results
brackenPreResults<-brackenPreResults[order(brackenPreResults$brackenPreTtest_p), , drop = F]
brackenPostResults<-brackenPostResults[order(brackenPostResults$brackenPostTtest_p), , drop = F]
amrPreResults<-amrPreResults[order(amrPreResults$amrPreTtest_p), , drop = F]
amrPostResults<-amrPostResults[order(amrPostResults$amrPostTtest_p), , drop = F]
rgiPreResults<-rgiPreResults[order(rgiPreResults$rgiPreTtest_p), , drop = F]
rgiPostResults<-rgiPostResults[order(rgiPostResults$rgiPostTtest_p), , drop = F]
vsearchPreResults<-vsearchPreResults[order(vsearchPreResults$vsearchPreTtest_p), , drop = F]
vsearchPostResults<-vsearchPostResults[order(vsearchPostResults$vsearchPostTtest_p), , drop = F]
genusPreResults<-genusPreResults[order(genusPreResults$genusPreTtest_p), , drop = F]
genusPostResults<-genusPostResults[order(genusPostResults$genusPostTtest_p), , drop = F]

#counts tables
brackenT<-read.csv("CountsTables/brackenFiltered.csv", header = T, row.names = 1, check.names = F)
amrT<-read.csv("CountsTables/amrFiltered.csv", header = T, row.names = 1, check.names = F)
rgiT<-read.csv("CountsTables/rgiFiltered.csv", header = T, row.names = 1, check.names = F)
vsearchT<-read.csv("CountsTables/vsearchFiltered.csv", header = T, row.names = 1, check.names = F)
genusT<-read.csv("CountsTables/genusFiltered.csv", header = T, row.names = 1, check.names = F)

patients<-unique(metaData$ID)

pres<-vector()
posts<-vector()
for (i in 1:length(patients)) {
  patient<-metaData[metaData$ID == patients[i], , drop = F]
  pres[i]<-rownames(patient)[patient$bins == "PRE"]
  posts[i]<-rownames(patient)[which.min(abs(patient$Timepoint - 28))]
}

brackenPre<-brackenT[, pres, drop = F]
brackenPost<-brackenT[, posts, drop = F]
amrPre<-amrT[, pres, drop = F]
amrPost<-amrT[, intersect(colnames(amrT), posts), drop = F]
rgiPre<-rgiT[, pres, drop = F]
rgiPost<-rgiT[, posts, drop = F]
vsearchPre<-vsearchT[, pres, drop = F]
vsearchPost<-vsearchT[, posts, drop = F]
genusPre<-genusT[, pres, drop = F]
genusPost<-genusT[, posts, drop = F]

#sort counts tables based on ttest results
brackenPre<-brackenPre[rownames(brackenPreResults), , drop = F]
brackenPost<-brackenPost[rownames(brackenPostResults), , drop = F]
amrPre<-amrPre[rownames(amrPreResults), , drop = F]
amrPost<-amrPost[rownames(amrPostResults), , drop = F]
rgiPre<-rgiPre[rownames(rgiPreResults), , drop = F]
rgiPost<-rgiPost[rownames(rgiPostResults), , drop = F]
vsearchPre<-vsearchPre[rownames(vsearchPreResults), , drop = F]
vsearchPost<-vsearchPost[rownames(vsearchPostResults), , drop = F]
genusPre<-genusPre[rownames(genusPreResults), , drop = F]
genusPost<-genusPost[rownames(genusPostResults), , drop = F]

metaBRACKENpre<-metaData[colnames(brackenPre), , drop = F]
metaBRACKENpost<-metaData[colnames(brackenPost), , drop = F]
metaAMRpre<-metaData[colnames(amrPre), , drop = F]
metaAMRpost<-metaData[colnames(amrPost), , drop = F]
metaRGIpre<-metaData[colnames(rgiPre), , drop = F]
metaRGIpost<-metaData[colnames(rgiPost), , drop = F]
metaVSEARCHpre<-metaData[colnames(vsearchPre), , drop = F]
metaVSEARCHpost<-metaData[colnames(vsearchPost), , drop = F]
metaGENUSpre<-metaData[colnames(genusPre), , drop = F]
metaGENUSpost<-metaData[colnames(genusPost), , drop = F]

#bracken pre
pdf("Plots/BrackenSpeciesPreInOutPatientSortedBoxPlots(Ttest).pdf", width=14, height=7)
par(mfrow=c(2, 4))
par(mar=c(5,6,4,1)+.1)
label = c("In", "Out")
for (i in 1:nrow(brackenT)){
  boxplot(unlist(brackenPre[i,])~metaBRACKENpre$ptInOut, outline = F,
          ylab=rownames(brackenPre)[i], xlab="Pre Patient", 
          names = label,
          cex.lab = 1.5, cex.axis = 1, cex.main = 1.8, 
          col=c("cornflowerblue", "darkorange"))
  stripchart(unlist(brackenPre[i,])~metaBRACKENpre$ptInOut, method="jitter", 
             vertical=T, pch=19, add=T)
  
  mtext(paste("Ttest P-value:",signif(brackenPreResults[rownames(brackenPre[i,]),1],digits = 3)),3, 0.9, cex = 0.6, adj = 0)
  mtext(paste("Adj. Ttest P-value:",signif(brackenPreResults[rownames(brackenPre[i,]),2],digits = 3)),3, 0.9, cex = 0.6, adj = 1)
}
dev.off()
#bracken post
pdf("Plots/BrackenSpeciesPostInOutPatientSortedBoxPlots(Ttest).pdf", width=14, height=7)
par(mfrow=c(2, 4))
par(mar=c(5,6,4,1)+.1)
label = c("In", "Out")
for (i in 1:nrow(brackenT)){
  boxplot(unlist(brackenPost[i,])~metaBRACKENpost$ptInOut, outline = F,
          ylab=rownames(brackenPost)[i], xlab="Post Patient", 
          names = label,
          cex.lab = 1.5, cex.axis = 1, cex.main = 1.8, 
          col=c("cornflowerblue", "darkorange"))
  stripchart(unlist(brackenPost[i,])~metaBRACKENpost$ptInOut, method="jitter", 
             vertical=T, pch=19, add=T)
  
  mtext(paste("Ttest P-value:",signif(brackenPostResults[rownames(brackenPost[i,]),1],digits = 3)),3, 0.9, cex = 0.6, adj = 0)
  mtext(paste("Adj. Ttest P-value:",signif(brackenPostResults[rownames(brackenPost[i,]),2],digits = 3)),3, 0.9, cex = 0.6, adj = 1)
}
dev.off()

#genus pre
pdf("Plots/BrackenGenusPreInOutPatientSortedBoxPlots(Ttest).pdf", width=14, height=7)
par(mfrow=c(2, 4))
par(mar=c(5,6,4,1)+.1)
label = c("In", "Out")
for (i in 1:nrow(genusT)){
  boxplot(unlist(genusPre[i,])~metaGENUSpre$ptInOut, outline = F,
          ylab=rownames(genusPre)[i], xlab="Pre Patient", 
          names = label,
          cex.lab = 1.5, cex.axis = 1, cex.main = 1.8, 
          col=c("cornflowerblue", "darkorange"))
  stripchart(unlist(genusPre[i,])~metaGENUSpre$ptInOut, method="jitter", 
             vertical=T, pch=19, add=T)
  
  mtext(paste("Ttest P-value:",signif(genusPreResults[rownames(genusPre[i,]),1],digits = 3)),3, 0.9, cex = 0.6, adj = 0)
  mtext(paste("Adj. Ttest P-value:",signif(genusPreResults[rownames(genusPre[i,]),2],digits = 3)),3, 0.9, cex = 0.6, adj = 1)
}
dev.off()
#genus post
pdf("Plots/BrackenGenusPostInOutPatientSortedBoxPlots(Ttest).pdf", width=14, height=7)
par(mfrow=c(2, 4))
par(mar=c(5,6,4,1)+.1)
label = c("In", "Out")
for (i in 1:nrow(genusT)){
  boxplot(unlist(genusPost[i,])~metaBRACKENpost$ptInOut, outline = F,
          ylab=rownames(genusPost)[i], xlab="Post Patient", 
          names = label,
          cex.lab = 1.5, cex.axis = 1, cex.main = 1.8, 
          col=c("cornflowerblue", "darkorange"))
  stripchart(unlist(genusPost[i,])~metaBRACKENpost$ptInOut, method="jitter", 
             vertical=T, pch=19, add=T)
  
  mtext(paste("Ttest P-value:",signif(genusPostResults[rownames(genusPost[i,]),1],digits = 3)),3, 0.9, cex = 0.6, adj = 0)
  mtext(paste("Adj. Ttest P-value:",signif(genusPostResults[rownames(genusPost[i,]),2],digits = 3)),3, 0.9, cex = 0.6, adj = 1)
}
dev.off()

#amr pre
pdf("Plots/AMRPreInOutPatientSortedBoxPlots(Ttest).pdf", width=14, height=7)
par(mfrow=c(2, 4))
par(mar=c(5,6,4,1)+.1)
label = c("In", "Out")
for (i in 1:nrow(amrT)){
  boxplot(unlist(amrPre[i,])~metaAMRpre$ptInOut, outline = F,
          ylab=rownames(amrPre)[i], xlab="Pre Patient", 
          names = label,
          cex.lab = 1.5, cex.axis = 1, cex.main = 1.8, 
          col=c("cornflowerblue", "darkorange"))
  stripchart(unlist(amrPre[i,])~metaAMRpre$ptInOut, method="jitter", 
             vertical=T, pch=19, add=T)
  
  mtext(paste("Ttest P-value:",signif(amrPreResults[rownames(amrPre[i,]),1],digits = 3)),3, 0.9, cex = 0.6, adj = 0)
  mtext(paste("Adj. Ttest P-value:",signif(amrPreResults[rownames(amrPre[i,]),2],digits = 3)),3, 0.9, cex = 0.6, adj = 1)
}
dev.off()
#amr post
pdf("Plots/AMRPostInOutPatientSortedBoxPlots(Ttest).pdf", width=14, height=7)
par(mfrow=c(2, 4))
par(mar=c(5,6,4,1)+.1)
label = c("In", "Out")
for (i in 1:nrow(amrT)){
  boxplot(unlist(amrPost[i,])~metaAMRpost$ptInOut, outline = F,
          ylab=rownames(amrPost)[i], xlab="Post Patient", 
          names = label,
          cex.lab = 1.5, cex.axis = 1, cex.main = 1.8, 
          col=c("cornflowerblue", "darkorange"))
  stripchart(unlist(amrPost[i,])~metaAMRpost$ptInOut, method="jitter", 
             vertical=T, pch=19, add=T)
  
  mtext(paste("Ttest P-value:",signif(amrPostResults[rownames(amrPost[i,]),1],digits = 3)),3, 0.9, cex = 0.6, adj = 0)
  mtext(paste("Adj. Ttest P-value:",signif(amrPostResults[rownames(amrPost[i,]),2],digits = 3)),3, 0.9, cex = 0.6, adj = 1)
}
dev.off()

#rgi pre
pdf("Plots/RGIPreInOutPatientSortedBoxPlots(Ttest).pdf", width=14, height=7)
par(mfrow=c(2, 4))
par(mar=c(5,6,4,1)+.1)
label = c("In", "Out")
for (i in 1:nrow(rgiT)){
  boxplot(unlist(rgiPre[i,])~metaRGIpre$ptInOut, outline = F,
          ylab=rownames(rgiPre)[i], xlab="Pre Patient", 
          names = label,
          cex.lab = 1.5, cex.axis = 1, cex.main = 1.8, 
          col=c("cornflowerblue", "darkorange"))
  stripchart(unlist(rgiPre[i,])~metaRGIpre$ptInOut, method="jitter", 
             vertical=T, pch=19, add=T)
  
  mtext(paste("Ttest P-value:",signif(rgiPreResults[rownames(rgiPre[i,]),1],digits = 3)),3, 0.9, cex = 0.6, adj = 0)
  mtext(paste("Adj. Ttest P-value:",signif(rgiPreResults[rownames(rgiPre[i,]),2],digits = 3)),3, 0.9, cex = 0.6, adj = 1)
}
dev.off()
#rgi post
pdf("Plots/RGIPostInOutPatientSortedBoxPlots(Ttest).pdf", width=14, height=7)
par(mfrow=c(2, 4))
par(mar=c(5,6,4,1)+.1)
label = c("In", "Out")
for (i in 1:nrow(rgiT)){
  boxplot(unlist(rgiPost[i,])~metaRGIpost$ptInOut, outline = F,
          ylab=rownames(rgiPost)[i], xlab="Post Patient", 
          names = label,
          cex.lab = 1.5, cex.axis = 1, cex.main = 1.8, 
          col=c("cornflowerblue", "darkorange"))
  stripchart(unlist(rgiPost[i,])~metaRGIpost$ptInOut, method="jitter", 
             vertical=T, pch=19, add=T)
  
  mtext(paste("Ttest P-value:",signif(rgiPostResults[rownames(rgiPost[i,]),1],digits = 3)),3, 0.9, cex = 0.6, adj = 0)
  mtext(paste("Adj. Ttest P-value:",signif(rgiPostResults[rownames(rgiPost[i,]),2],digits = 3)),3, 0.9, cex = 0.6, adj = 1)
}
dev.off()

#vsearch pre
pdf("Plots/VsearchPreInOutPatientSortedBoxPlots(Ttest).pdf", width=14, height=7)
par(mfrow=c(2, 4))
par(mar=c(5,6,4,1)+.1)
label = c("In", "Out")
for (i in 1:nrow(vsearchT)){
  boxplot(unlist(vsearchPre[i,])~metaVSEARCHpre$ptInOut, outline = F,
          ylab=rownames(vsearchPre)[i], xlab="Pre Patient", 
          names = label,
          cex.lab = 1.5, cex.axis = 1, cex.main = 1.8, 
          col=c("cornflowerblue", "darkorange"))
  stripchart(unlist(vsearchPre[i,])~metaVSEARCHpre$ptInOut, method="jitter", 
             vertical=T, pch=19, add=T)
  
  mtext(paste("Ttest P-value:",signif(vsearchPreResults[rownames(vsearchPre[i,]),1],digits = 3)),3, 0.9, cex = 0.6, adj = 0)
  mtext(paste("Adj. Ttest P-value:",signif(vsearchPreResults[rownames(vsearchPre[i,]),2],digits = 3)),3, 0.9, cex = 0.6, adj = 1)
}
dev.off()
#vsearch post
pdf("Plots/VsearchPostInOutPatientSortedBoxPlots(Ttest).pdf", width=14, height=7)
par(mfrow=c(2, 4))
par(mar=c(5,6,4,1)+.1)
label = c("In", "Out")
for (i in 1:nrow(vsearchT)){
  boxplot(unlist(vsearchPost[i,])~metaVSEARCHpost$ptInOut, outline = F,
          ylab=rownames(vsearchPost)[i], xlab="Post Patient", 
          names = label,
          cex.lab = 1.5, cex.axis = 1, cex.main = 1.8, 
          col=c("cornflowerblue", "darkorange"))
  stripchart(unlist(vsearchPost[i,])~metaVSEARCHpost$ptInOut, method="jitter", 
             vertical=T, pch=19, add=T)
  
  mtext(paste("Ttest P-value:",signif(vsearchPostResults[rownames(vsearchPost[i,]),1],digits = 3)),3, 0.9, cex = 0.6, adj = 0)
  mtext(paste("Adj. Ttest P-value:",signif(vsearchPostResults[rownames(vsearchPost[i,]),2],digits = 3)),3, 0.9, cex = 0.6, adj = 1)
}
dev.off()