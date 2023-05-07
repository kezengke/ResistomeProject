#total counts for each patient scatter plots
rm(list = ls())
library(stringr)

metaData<-read.csv("metaWithBins.csv", header = T, row.names = 1)

amrT<-read.csv("CountsTables/amrProcessed.csv", header = T, row.names = 1, check.names = F)
metaAMR<-metaData[colnames(amrT), , drop = F]
ID<-unique(metaAMR$PID)

#AMR
pdf("Plots/AMRTotalCountsForEachSample.pdf", width=10, height=15)
par(mfrow=c(3,2))
for (i in 1:length(ID)) {
  meta<-metaAMR[metaAMR$PID == ID[i], ,drop = F]
  samples<-rownames(meta)
  counts<-amrT[, samples, drop = F]
  totalCounts<-colSums(counts)
  
  plot(meta$Timepoint, totalCounts, col = "coral3", pch = 19, 
       xlab = "Days", ylab = "AMR Counts",
       main = paste0("D", ID[i]))
  
}
dev.off()

#RGI
rgiT<-read.csv("CountsTables/rgiProcessed.csv", header = T, row.names = 1, check.names = F)
metaRGI<-metaData[colnames(rgiT), , drop = F]
ID<-unique(metaRGI$PID)

pdf("Plots/RGITotalCountsForEachSample.pdf", width=10, height=15)
par(mfrow=c(3,2))
for (i in 1:length(ID)) {
  meta<-metaRGI[metaRGI$PID == ID[i], ,drop = F]
  samples<-rownames(meta)
  counts<-rgiT[, samples, drop = F]
  totalCounts<-colSums(counts)
  
  plot(meta$Timepoint, totalCounts, col = "cornflowerblue", pch = 19, 
       xlab = "Days", ylab = "RGI Counts",
       main = paste0("D", ID[i]))
  
}
dev.off()

#vsearch
vsearchT<-read.csv("CountsTables/vsearchProcessed.csv", header = T, row.names = 1, check.names = F)
metaVSEARCH<-metaData[colnames(vsearchT), , drop = F]
ID<-unique(metaVSEARCH$PID)

pdf("Plots/vsearchTotalCountsForEachSample.pdf", width=10, height=15)
par(mfrow=c(3,2))
for (i in 1:length(ID)) {
  meta<-metaVSEARCH[metaVSEARCH$PID == ID[i], ,drop = F]
  samples<-rownames(meta)
  counts<-vsearchT[, samples, drop = F]
  totalCounts<-colSums(counts)
  
  plot(meta$Timepoint, totalCounts, col = "olivedrab4", pch = 19, 
       xlab = "Days", ylab = "vsearch Counts",
       main = paste0("D", ID[i]))
  
}
dev.off()