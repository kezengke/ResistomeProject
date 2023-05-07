#days vs gene type richness for each patient scatter plots
rm(list = ls())
library(stringr)

metaData<-read.csv("metaWithBins.csv", header = T, row.names = 1)

#AMR
amrT<-read.csv("CountsTables/amrProcessed.csv", header = T, row.names = 1, check.names = F)
metaAMR<-metaData[colnames(amrT), , drop = F]
ID<-unique(metaAMR$PID)

pdf("Plots/AMRGeneTypeRichnessForEachSample.pdf", width=10, height=15)
par(mfrow=c(3,2))
for (i in 1:length(ID)) {
  meta<-metaAMR[metaAMR$PID == ID[i], ,drop = F]
  samples<-rownames(meta)
  counts<-amrT[, samples, drop = F]
  geneTypes<-colSums(counts>0)
  
  plot(meta$Timepoint, geneTypes, col = "coral3", pch = 19, 
       xlab = "Days", ylab = "Types of Genes",
       main = paste0("D", ID[i], " richness(AMR)"))
}
dev.off()

#RGI
rgiT<-read.csv("CountsTables/rgiProcessed.csv", header = T, row.names = 1, check.names = F)
metaRGI<-metaData[colnames(rgiT), , drop = F]
ID<-unique(metaRGI$PID)

pdf("Plots/RGIGeneTypeRichnessForEachSample.pdf", width=10, height=15)
par(mfrow=c(3,2))
for (i in 1:length(ID)) {
  meta<-metaRGI[metaRGI$PID == ID[i], ,drop = F]
  samples<-rownames(meta)
  counts<-rgiT[, samples, drop = F]
  geneTypes<-colSums(counts>0)
  
  plot(meta$Timepoint, geneTypes, col = "cornflowerblue", pch = 19, 
       xlab = "Days", ylab = "Types of Genes",
       main = paste0("D", ID[i], " richness(RGI)"))
}
dev.off()

#VSEARCH
vsearchT<-read.csv("CountsTables/vsearchProcessed.csv", header = T, row.names = 1, check.names = F)
metaVSEARCH<-metaData[colnames(vsearchT), , drop = F]
ID<-unique(metaVSEARCH$PID)

pdf("Plots/VsearchGeneTypeRichnessForEachSample.pdf", width=10, height=15)
par(mfrow=c(3,2))
for (i in 1:length(ID)) {
  meta<-metaVSEARCH[metaVSEARCH$PID == ID[i], ,drop = F]
  samples<-rownames(meta)
  counts<-vsearchT[, samples, drop = F]
  geneTypes<-colSums(counts>0)
  
  plot(meta$Timepoint, geneTypes, col = "olivedrab4", pch = 19, 
       xlab = "Days", ylab = "Types of Genes",
       main = paste0("D", ID[i], " richness(VSEARCH)"))
}
dev.off()
