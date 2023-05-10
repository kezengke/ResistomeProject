rm(list = ls())
library("vegan")
library("RColorBrewer")

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

circleCol<-brewer.pal(length(unique(metaBRACKEN$bins)), "Paired")
cols<-circleCol[factor(metaBRACKEN$bins, levels = c("PRE", "D0", "D7", "D14", "D21", "D28", "D35", "D60", "D100", "D180", "D365", "D730"))]

IDtypes<-unique(metaBRACKEN$ID)
pdf("Plots/MDSPlotsForEachPatient(bracken).pdf", width=12, height=18)
par(mfrow=c(3,2))
par(mar=c(5,6,4,1)+.1)
for (i in 1:length(IDtypes)) {
  mdsT<-brackenT[, colnames(brackenT)[metaBRACKEN$ID == IDtypes[i]]]
  if(ncol(mdsT)<3)
    next
  MDS<-capscale(t(mdsT)~1,distance = "bray")
  percentVariance<-MDS$CA$eig/sum(eigenvals(MDS))*100
  statusPlot<-ordiplot(MDS,choices = c(1,2),type="none",cex.lab = 1,
                       xlab = paste("MDS1  ", format(percentVariance[1], digits = 4), "%", sep = ""),
                       ylab = paste("MDS2  ", format(percentVariance[2], digits = 4), "%", sep = ""),
                       main = paste("Bracken-Species", IDtypes[i]),
                       xlim = c(-2.5, 1.5), ylim = c(-2, 2))
  points(statusPlot,"sites", pch = 19, cex = 2.5, col = adjustcolor(cols[which(metaBRACKEN$ID == IDtypes[i])], alpha.f = 0.5))
  legend("topright", c("PRE", "D0", "D7", "D14", "D21", "D28", "D35", "D60", "D100", "D180", "D365", "D730"), col = circleCol[1:12], cex = 1, pch = 16, bty = "n")
}
dev.off()

for (i in 1:length(IDtypes)) {
  mdsT<-brackenT[, colnames(brackenT)[metaBRACKEN$ID == IDtypes[i]]]
  if(ncol(mdsT)<3)
    next
  MDS<-capscale(t(mdsT)~1,distance = "bray")
  percentVariance<-MDS$CA$eig/sum(eigenvals(MDS))*100
  pdf(paste0("MDSPlotsForEachPatient/Bracken/Bracken-Species", IDtypes[i], ".pdf"), width=6, height=6)
  statusPlot<-ordiplot(MDS,choices = c(1,2),type="none",cex.lab = 1,
                       xlab = paste("MDS1  ", format(percentVariance[1], digits = 4), "%", sep = ""),
                       ylab = paste("MDS2  ", format(percentVariance[2], digits = 4), "%", sep = ""),
                       main = paste("Bracken-Species", IDtypes[i]),
                       xlim = c(-2.5, 1.5), ylim = c(-2, 2))
  points(statusPlot,"sites", pch = 19, cex = 2.5, col = adjustcolor(cols[which(metaBRACKEN$ID == IDtypes[i])], alpha.f = 0.5))
  legend("topright", c("PRE", "D0", "D7", "D14", "D21", "D28", "D35", "D60", "D100", "D180", "D365", "D730"), col = circleCol[1:12], cex = 1, pch = 16, bty = "n")
  dev.off()
}

circleCol<-brewer.pal(length(unique(metaAMR$bins)), "Paired")
cols<-circleCol[factor(metaAMR$bins, levels = c("PRE", "D0", "D7", "D14", "D21", "D28", "D35", "D60", "D100", "D180", "D365", "D730"))]

IDtypes<-unique(metaAMR$ID)
pdf("Plots/MDSPlotsForEachPatient(AMR).pdf", width=12, height=18)
par(mfrow=c(3,2))
par(mar=c(5,6,4,1)+.1)
for (i in 1:length(IDtypes)) {
  mdsT<-amrT[, colnames(amrT)[metaAMR$ID == IDtypes[i]]]
  if(ncol(mdsT)<3)
    next
  MDS<-capscale(t(mdsT)~1,distance = "bray")
  percentVariance<-MDS$CA$eig/sum(eigenvals(MDS))*100
  statusPlot<-ordiplot(MDS,choices = c(1,2),type="none",cex.lab = 1,
                       xlab = paste("MDS1  ", format(percentVariance[1], digits = 4), "%", sep = ""),
                       ylab = paste("MDS2  ", format(percentVariance[2], digits = 4), "%", sep = ""),
                       main = paste("AMR", IDtypes[i]),
                       xlim = c(-3, 1.5), ylim = c(-2, 2.5))
  points(statusPlot,"sites", pch = 19, cex = 2.5, col = adjustcolor(cols[which(metaAMR$ID == IDtypes[i])], alpha.f = 0.5))
  legend("topright", c("PRE", "D0", "D7", "D14", "D21", "D28", "D35", "D60", "D100", "D180", "D365", "D730"), col = circleCol[1:12], cex = 1, pch = 16, bty = "n")
}
dev.off()

for (i in 1:length(IDtypes)) {
  mdsT<-amrT[, colnames(amrT)[metaAMR$ID == IDtypes[i]]]
  if(ncol(mdsT)<3)
    next
  MDS<-capscale(t(mdsT)~1,distance = "bray")
  percentVariance<-MDS$CA$eig/sum(eigenvals(MDS))*100
  pdf(paste0("MDSPlotsForEachPatient/AMR/AMR", IDtypes[i], ".pdf"), width=6, height=6)
  statusPlot<-ordiplot(MDS,choices = c(1,2),type="none",cex.lab = 1,
                       xlab = paste("MDS1  ", format(percentVariance[1], digits = 4), "%", sep = ""),
                       ylab = paste("MDS2  ", format(percentVariance[2], digits = 4), "%", sep = ""),
                       main = paste("AMR", IDtypes[i]),
                       xlim = c(-3, 1.5), ylim = c(-2, 2.5))
  points(statusPlot,"sites", pch = 19, cex = 2.5, col = adjustcolor(cols[which(metaAMR$ID == IDtypes[i])], alpha.f = 0.5))
  legend("topright", c("PRE", "D0", "D7", "D14", "D21", "D28", "D35", "D60", "D100", "D180", "D365", "D730"), col = circleCol[1:12], cex = 1, pch = 16, bty = "n")
  dev.off()
}

circleCol<-brewer.pal(length(unique(metaRGI$bins)), "Paired")
cols<-circleCol[factor(metaRGI$bins, levels = c("PRE", "D0", "D7", "D14", "D21", "D28", "D35", "D60", "D100", "D180", "D365", "D730"))]

IDtypes<-unique(metaRGI$ID)
pdf("Plots/MDSPlotsForEachPatient(RGI).pdf", width=12, height=18)
par(mfrow=c(3,2))
par(mar=c(5,6,4,1)+.1)
for (i in 1:length(IDtypes)) {
  mdsT<-rgiT[, colnames(rgiT)[metaRGI$ID == IDtypes[i]]]
  if(ncol(mdsT)<3)
    next
  MDS<-capscale(t(mdsT)~1,distance = "bray")
  percentVariance<-MDS$CA$eig/sum(eigenvals(MDS))*100
  statusPlot<-ordiplot(MDS,choices = c(1,2),type="none",cex.lab = 1,
                       xlab = paste("MDS1  ", format(percentVariance[1], digits = 4), "%", sep = ""),
                       ylab = paste("MDS2  ", format(percentVariance[2], digits = 4), "%", sep = ""),
                       main = paste("RGI", IDtypes[i]),
                       xlim = c(-2, 2), ylim = c(-2, 2))
  points(statusPlot,"sites", pch = 19, cex = 2.5, col = adjustcolor(cols[which(metaRGI$ID == IDtypes[i])], alpha.f = 0.5))
  legend("topright", c("PRE", "D0", "D7", "D14", "D21", "D28", "D35", "D60", "D100", "D180", "D365", "D730"), col = circleCol[1:13], cex = 1, pch = 16, bty = "n")
}
dev.off()

for (i in 1:length(IDtypes)) {
  mdsT<-rgiT[, colnames(rgiT)[metaRGI$ID == IDtypes[i]]]
  if(ncol(mdsT)<3)
    next
  MDS<-capscale(t(mdsT)~1,distance = "bray")
  percentVariance<-MDS$CA$eig/sum(eigenvals(MDS))*100
  pdf(paste0("MDSPlotsForEachPatient/RGI/RGI", IDtypes[i], ".pdf"), width=6, height=6)
  statusPlot<-ordiplot(MDS,choices = c(1,2),type="none",cex.lab = 1,
                       xlab = paste("MDS1  ", format(percentVariance[1], digits = 4), "%", sep = ""),
                       ylab = paste("MDS2  ", format(percentVariance[2], digits = 4), "%", sep = ""),
                       main = paste("RGI", IDtypes[i]),
                       xlim = c(-2, 2), ylim = c(-2, 2))
  points(statusPlot,"sites", pch = 19, cex = 2.5, col = adjustcolor(cols[which(metaRGI$ID == IDtypes[i])], alpha.f = 0.5))
  legend("topright", c("PRE", "D0", "D7", "D14", "D21", "D28", "D35", "D60", "D100", "D180", "D365", "D730"), col = circleCol[1:12], cex = 1, pch = 16, bty = "n")
  dev.off()
}

circleCol<-brewer.pal(length(unique(metaVSEARCH$bins)), "Paired")
cols<-circleCol[factor(metaVSEARCH$bins, levels = c("PRE", "D0", "D7", "D14", "D21", "D28", "D35", "D60", "D100", "D180", "D365", "D730"))]

IDtypes<-unique(metaVSEARCH$ID)
pdf("Plots/MDSPlotsForEachPatient(vsearch).pdf", width=12, height=18)
par(mfrow=c(3,2))
par(mar=c(5,6,4,1)+.1)
for (i in 1:length(IDtypes)) {
  mdsT<-vsearchT[, colnames(vsearchT)[metaVSEARCH$ID == IDtypes[i]]]
  if(ncol(mdsT)<3)
    next
  MDS<-capscale(t(mdsT)~1,distance = "bray")
  percentVariance<-MDS$CA$eig/sum(eigenvals(MDS))*100
  statusPlot<-ordiplot(MDS,choices = c(1,2),type="none",cex.lab = 1,
                       xlab = paste("MDS1  ", format(percentVariance[1], digits = 4), "%", sep = ""),
                       ylab = paste("MDS2  ", format(percentVariance[2], digits = 4), "%", sep = ""),
                       main = paste("vsearch", IDtypes[i]),
                       xlim = c(-1.5, 2.5), ylim = c(-2, 2))
  points(statusPlot,"sites", pch = 19, cex = 2.5, col = adjustcolor(cols[which(metaVSEARCH$ID == IDtypes[i])], alpha.f = 0.5))
  legend("topright", c("PRE", "D0", "D7", "D14", "D21", "D28", "D35", "D60", "D100", "D180", "D365", "D730"), col = circleCol[1:12], cex = 1, pch = 16, bty = "n")
}
dev.off()

for (i in 1:length(IDtypes)) {
  mdsT<-vsearchT[, colnames(vsearchT)[metaVSEARCH$ID == IDtypes[i]]]
  if(ncol(mdsT)<3)
    next
  MDS<-capscale(t(mdsT)~1,distance = "bray")
  percentVariance<-MDS$CA$eig/sum(eigenvals(MDS))*100
  pdf(paste0("MDSPlotsForEachPatient/vsearch/Vsearch", IDtypes[i], ".pdf"), width=6, height=6)
  statusPlot<-ordiplot(MDS,choices = c(1,2),type="none",cex.lab = 1,
                       xlab = paste("MDS1  ", format(percentVariance[1], digits = 4), "%", sep = ""),
                       ylab = paste("MDS2  ", format(percentVariance[2], digits = 4), "%", sep = ""),
                       main = paste("vsearch", IDtypes[i]),
                       xlim = c(-1.5, 2.5), ylim = c(-2, 2))
  points(statusPlot,"sites", pch = 19, cex = 2.5, col = adjustcolor(cols[which(metaVSEARCH$ID == IDtypes[i])], alpha.f = 0.5))
  legend("topright", c("PRE", "D0", "D7", "D14", "D21", "D28", "D35", "D60", "D100", "D180", "D365", "D730"), col = circleCol[1:12], cex = 1, pch = 16, bty = "n")
  dev.off()
}