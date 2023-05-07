#correlation histograms for each patient total counts and days for each gene counts table
rm(list = ls())
library(stringr)
pdf("Plots/CorrelationHistograms.pdf", width=10, height=15)
par(mfrow=c(3,2))

metaData<-read.csv("metaWithBins.csv", header = T, row.names = 1)

#AMR
amrT<-read.csv("CountsTables/amrProcessed.csv", header = T, row.names = 1, check.names = F)
metaAMR<-metaData[colnames(amrT), , drop = F]
ID<-unique(metaAMR$PID)

pearsonCor<-c()
spearmanRho<-c()

for (i in 1:length(ID)) {
  meta<-metaAMR[metaAMR$PID == ID[i], ,drop = F]
  if(nrow(meta)<3)
    next
  samples<-rownames(meta)
  counts<-amrT[, samples, drop = F]
  totalCounts<-colSums(counts)
  
  pearsonCor<-c(pearsonCor, cor.test(meta$Timepoint, totalCounts, method = 'pearson')$estimate)
  spearmanRho<-c(spearmanRho, cor.test(meta$Timepoint, totalCounts, method = 'spearman')$estimate)
  
}

hist(pearsonCor, breaks=seq(-1, 1, 0.1), col = "coral3", 
     xlab = "Pearson's Cor", main="AMR", 
     cex.lab = 1.5, cex.main = 1.7, cex.axis = 1.4)
hist(spearmanRho, breaks=seq(-1, 1, 0.1), col = "coral3", 
     xlab = "Spearman's Rho", main="AMR", 
     cex.lab = 1.5, cex.main = 1.7, cex.axis = 1.4)

#RGI
rgiT<-read.csv("CountsTables/rgiProcessed.csv", header = T, row.names = 1, check.names = F)
metaRGI<-metaData[colnames(rgiT), , drop = F]
ID<-unique(metaRGI$PID)

pearsonCor<-c()
spearmanRho<-c()

for (i in 1:length(ID)) {
  meta<-metaRGI[metaRGI$PID == ID[i], ,drop = F]
  if(nrow(meta)<3)
    next
  samples<-rownames(meta)
  counts<-rgiT[, samples, drop = F]
  totalCounts<-colSums(counts)
  
  pearsonCor<-c(pearsonCor, cor.test(meta$Timepoint, totalCounts, method = 'pearson')$estimate)
  spearmanRho<-c(spearmanRho, cor.test(meta$Timepoint, totalCounts, method = 'spearman')$estimate)
  
}

hist(pearsonCor, breaks=seq(-1, 1, 0.1), col = "cornflowerblue", 
     xlab = "Pearson's Cor", main="RGI", 
     cex.lab = 1.5, cex.main = 1.7, cex.axis = 1.4)
hist(spearmanRho, breaks=seq(-1, 1, 0.1), col = "cornflowerblue", 
     xlab = "Spearman's Rho", main="RGI", 
     cex.lab = 1.5, cex.main = 1.7, cex.axis = 1.4)

#VSEARCH
vsearchT<-read.csv("CountsTables/vsearchProcessed.csv", header = T, row.names = 1, check.names = F)
metaVSEARCH<-metaData[colnames(vsearchT), , drop = F]
ID<-unique(metaVSEARCH$PID)

pearsonCor<-c()
spearmanRho<-c()

for (i in 1:length(ID)) {
  meta<-metaVSEARCH[metaVSEARCH$PID == ID[i], ,drop = F]
  if(nrow(meta)<3)
    next
  samples<-rownames(meta)
  counts<-vsearchT[, samples, drop = F]
  totalCounts<-colSums(counts)
  
  pearsonCor<-c(pearsonCor, cor.test(meta$Timepoint, totalCounts, method = 'pearson')$estimate)
  spearmanRho<-c(spearmanRho, cor.test(meta$Timepoint, totalCounts, method = 'spearman')$estimate)
  
}

hist(pearsonCor, breaks=seq(-1, 1, 0.1), col = "olivedrab4", 
     xlab = "Pearson's Cor", main="VSEARCH", 
     cex.lab = 1.5, cex.main = 1.7, cex.axis = 1.4)
hist(spearmanRho, breaks=seq(-1, 1, 0.1), col = "olivedrab4", 
     xlab = "Spearman's Rho", main="VSEARCH", 
     cex.lab = 1.5, cex.main = 1.7, cex.axis = 1.4)

dev.off()
