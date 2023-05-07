rm(list = ls())
library("vegan")

metaData<-read.csv("metaWithBins.csv", header = T, row.names = 1)

#counts tables
brackenT<-read.csv("CountsTables/brackenNormalized.csv", header = T, row.names = 1, check.names = F)
amrT<-read.csv("CountsTables/amrNormalized.csv", header = T, row.names = 1, check.names = F)
rgiT<-read.csv("CountsTables/rgiNormalized.csv", header = T, row.names = 1, check.names = F)
vsearchT<-read.csv("CountsTables/vsearchNormalized.csv", header = T, row.names = 1, check.names = F)

metaBRACKEN<-metaData[colnames(brackenT), , drop = F]
metaAMR<-metaData[colnames(amrT), , drop = F]
metaRGI<-metaData[colnames(rgiT), , drop = F]
metaVSEARCH<-metaData[colnames(vsearchT), , drop = F]

TimePointtypes<-unique(metaBRACKEN$bins)
for (i in 1:length(TimePointtypes)) {
  mdsT<-brackenT[, colnames(brackenT)[metaBRACKEN$bins == TimePointtypes[i]]]
  if(ncol(mdsT)<3)
    next
  MDS<-capscale(t(mdsT)~1,distance = "bray")
  percentVariance<-MDS$CA$eig/sum(eigenvals(MDS))*100
  pdf(paste("MDSPlotsForEachTimepoint/Bracken/Bracken-Species", TimePointtypes[i], ".pdf"), width=6, height=6)
  statusPlot<-ordiplot(MDS,choices = c(1,2),type="none",cex.lab = 1,
                       xlab = paste("MDS1  ", format(percentVariance[1], digits = 4), "%", sep = ""),
                       ylab = paste("MDS2  ", format(percentVariance[2], digits = 4), "%", sep = ""),
                       main = paste0("Bracken-Species ", TimePointtypes[i], " (n=", ncol(mdsT), ")"),
                       xlim = c(-2.5, 2.5), ylim = c(-2, 2))
  points(statusPlot,"sites", pch = 19, cex = 2.5, col = adjustcolor("tan2", alpha.f = 0.5))
  text(statusPlot, "sites", labels = metaBRACKEN[colnames(mdsT), 5], cex = 0.6, pos = 4, col = "tan2")
  dev.off()
}

TimePointtypes<-unique(metaAMR$bins)
for (i in 1:length(TimePointtypes)) {
  mdsT<-amrT[, colnames(amrT)[metaAMR$bins == TimePointtypes[i]]]
  if(ncol(mdsT)<3)
    next
  MDS<-capscale(t(mdsT)~1,distance = "bray")
  percentVariance<-MDS$CA$eig/sum(eigenvals(MDS))*100
  pdf(paste("MDSPlotsForEachTimepoint/AMR/AMR", TimePointtypes[i], ".pdf"), width=6, height=6)
  statusPlot<-ordiplot(MDS,choices = c(1,2),type="none",cex.lab = 1,
                       xlab = paste("MDS1  ", format(percentVariance[1], digits = 4), "%", sep = ""),
                       ylab = paste("MDS2  ", format(percentVariance[2], digits = 4), "%", sep = ""),
                       main = paste0("AMR ", TimePointtypes[i], " (n=", ncol(mdsT), ")"),
                       xlim = c(-3, 2.5), ylim = c(-2, 2.5))
  points(statusPlot,"sites", pch = 19, cex = 2.5, col = adjustcolor("coral3", alpha.f = 0.5))
  text(statusPlot, "sites", labels = metaAMR[colnames(mdsT), 5], cex = 0.6, pos = 4, col = "coral3")
  dev.off()
}

TimePointtypes<-unique(metaRGI$bins)
for (i in 1:length(TimePointtypes)) {
  mdsT<-rgiT[, colnames(rgiT)[metaRGI$bins == TimePointtypes[i]]]
  if(ncol(mdsT)<3)
    next
  MDS<-capscale(t(mdsT)~1,distance = "bray")
  percentVariance<-MDS$CA$eig/sum(eigenvals(MDS))*100
  pdf(paste("MDSPlotsForEachTimepoint/RGI/RGI", TimePointtypes[i], ".pdf"), width=6, height=6)
  statusPlot<-ordiplot(MDS,choices = c(1,2),type="none",cex.lab = 1,
                       xlab = paste("MDS1  ", format(percentVariance[1], digits = 4), "%", sep = ""),
                       ylab = paste("MDS2  ", format(percentVariance[2], digits = 4), "%", sep = ""),
                       main = paste0("RGI ", TimePointtypes[i], " (n=", ncol(mdsT), ")"),
                       xlim = c(-2, 2), ylim = c(-2, 3))
  points(statusPlot,"sites", pch = 19, cex = 2.5, col = adjustcolor("cornflowerblue", alpha.f = 0.5))
  text(statusPlot, "sites", labels = metaRGI[colnames(mdsT), 5], cex = 0.6, pos = 4, col = "cornflowerblue")
  dev.off()
}

TimePointtypes<-unique(metaVSEARCH$bins)
for (i in 1:length(TimePointtypes)) {
  mdsT<-vsearchT[, colnames(vsearchT)[metaVSEARCH$bins == TimePointtypes[i]]]
  if(ncol(mdsT)<3)
    next
  MDS<-capscale(t(mdsT)~1,distance = "bray")
  percentVariance<-MDS$CA$eig/sum(eigenvals(MDS))*100
  pdf(paste("MDSPlotsForEachTimepoint/vsearch/vsearch", TimePointtypes[i], ".pdf"), width=6, height=6)
  statusPlot<-ordiplot(MDS,choices = c(1,2),type="none",cex.lab = 1,
                       xlab = paste("MDS1  ", format(percentVariance[1], digits = 4), "%", sep = ""),
                       ylab = paste("MDS2  ", format(percentVariance[2], digits = 4), "%", sep = ""),
                       main = paste0("vsearch ", TimePointtypes[i], " (n=", ncol(mdsT), ")"),
                       xlim = c(-1.5, 2.5), ylim = c(-2, 2))
  points(statusPlot,"sites", pch = 19, cex = 2.5, col = adjustcolor("olivedrab4", alpha.f = 0.5))
  text(statusPlot, "sites", labels = metaVSEARCH[colnames(mdsT), 5], cex = 0.6, pos = 4, col = "olivedrab4")
  dev.off()
}
