#MDS plots for in and out patients for every time frame.
rm(list = ls())
library("vegan")

metaData<-read.csv("metaWithBins.csv", header = T, row.names = 1)

#counts tables
brackenT<-read.csv("CountsTables/brackenFiltered.csv", header = T, row.names = 1, check.names = F)
amrT<-read.csv("CountsTables/amrFiltered.csv", header = T, row.names = 1, check.names = F)
rgiT<-read.csv("CountsTables/rgiFiltered.csv", header = T, row.names = 1, check.names = F)
vsearchT<-read.csv("CountsTables/vsearchFiltered.csv", header = T, row.names = 1, check.names = F)

timeFrames<-unique(metaData$bins)

pdf("Plots/MDSForInAndOutPatients(AllTimeFrames).pdf", width=24, height=30)
par(mfrow=c(5,4))
par(mar=c(5,6,4,1)+.1)
for (i in 1:length(timeFrames)) {
  samples<-rownames(metaData)[metaData$bins == timeFrames[i]]
  bracken<-brackenT[, samples, drop = F]
  amr<-amrT[, intersect(colnames(amrT), samples), drop = F]
  rgi<-rgiT[, samples, drop = F]
  vsearch<-vsearchT[, samples, drop = F]
  
  metaBRACKEN<-metaData[colnames(bracken), , drop = F]
  metaAMR<-metaData[colnames(amr), , drop = F]
  metaRGI<-metaData[colnames(rgi), , drop = F]
  metaVSEARCH<-metaData[colnames(vsearch), , drop = F]
  
  circleCol<-c("cornflowerblue", "darkorange")
  #Bracken
  cols<-circleCol[factor(metaBRACKEN$ptInOut, levels = c("Outpatient", "Inpatient"))]
  
  MDS<-capscale(t(bracken)~1,distance = "bray")
  percentVariance<-MDS$CA$eig/sum(eigenvals(MDS))*100
  pval<-adonis2(t(bracken)~metaBRACKEN$ptInOut, method="bray")$`Pr(>F)`[1]
  statusPlot<-ordiplot(MDS,choices = c(1,2),type="none",cex.lab = 1,
                       xlab = paste("MDS1  ", format(percentVariance[1], digits = 4), "%", sep = ""),
                       ylab = paste("MDS2  ", format(percentVariance[2], digits = 4), "%", sep = ""),
                       main = paste(timeFrames[i], "Bracken-Species n =", ncol(bracken), "\nP-value:",pval))
  points(statusPlot,"sites", pch = 19, cex = 2.5, col = adjustcolor(cols, alpha.f = 0.5))
  ordiellipse(statusPlot, metaBRACKEN$ptInOut, kind="se", conf=0.95, lwd=4, draw = "lines", col=circleCol[1],show.groups="Outpatient",label=T,font=2,cex=1) 
  ordiellipse(statusPlot, metaBRACKEN$ptInOut, kind="se", conf=0.95, lwd=4, draw = "lines", col=circleCol[2],show.groups="Inpatient",label=T,font=2,cex=1) 
  legend("topright", c("Outpatient", "Inpatient"), col = circleCol[1:2], cex = 1, pch = 16, bty = "n")
  
  #AMR
  cols<-circleCol[factor(metaAMR$ptInOut, levels = c("Outpatient", "Inpatient"))]
  
  MDS<-capscale(t(amr)~1,distance = "bray")
  percentVariance<-MDS$CA$eig/sum(eigenvals(MDS))*100
  pval<-adonis2(t(amr)~metaAMR$ptInOut, method="bray")$`Pr(>F)`[1]
  statusPlot<-ordiplot(MDS,choices=c(1,2),type="none",cex.lab=1,
                       xlab=paste("MDS1  ", format(percentVariance[1], digits = 4), "%", sep = ""),
                       ylab=paste("MDS2  ", format(percentVariance[2], digits = 4), "%", sep = ""),
                       main=paste(timeFrames[i], "AMR n =", ncol(amr), "\nP-value:",pval))
  points(statusPlot,"sites", pch=19, cex=2.5, col=adjustcolor(cols, alpha.f = 0.5))
  ordiellipse(statusPlot, metaAMR$ptInOut, kind="se", conf=0.95, lwd=4, draw = "lines", col=circleCol[1],show.groups="Outpatient",label=T,font=2,cex=1) 
  ordiellipse(statusPlot, metaAMR$ptInOut, kind="se", conf=0.95, lwd=4, draw = "lines", col=circleCol[2],show.groups="Inpatient",label=T,font=2,cex=1) 
  legend("topright", c("Outpatient", "Inpatient"), col = circleCol[1:2], cex = 1, pch = 16, bty = "n")
  
  #RGI
  cols<-circleCol[factor(metaRGI$ptInOut, levels = c("Outpatient", "Inpatient"))]
  
  MDS<-capscale(t(rgi)~1,distance = "bray")
  percentVariance<-MDS$CA$eig/sum(eigenvals(MDS))*100
  pval<-adonis2(t(rgi)~metaRGI$ptInOut, method="bray")$`Pr(>F)`[1]
  statusPlot<-ordiplot(MDS,choices=c(1,2),type="none",cex.lab=1,
                       xlab=paste("MDS1  ", format(percentVariance[1], digits = 4), "%", sep = ""),
                       ylab=paste("MDS2  ", format(percentVariance[2], digits = 4), "%", sep = ""),
                       main=paste(timeFrames[i], "RGI n =", ncol(rgi), "\nP-value:",pval))
  points(statusPlot,"sites", pch=19, cex=2.5, col=adjustcolor(cols, alpha.f = 0.5))
  ordiellipse(statusPlot, metaRGI$ptInOut, kind="se", conf=0.95, lwd=4, draw = "lines", col=circleCol[1],show.groups="Outpatient",label=T,font=2,cex=1) 
  ordiellipse(statusPlot, metaRGI$ptInOut, kind="se", conf=0.95, lwd=4, draw = "lines", col=circleCol[2],show.groups="Inpatient",label=T,font=2,cex=1) 
  legend("topright", c("Outpatient", "Inpatient"), col = circleCol[1:2], cex = 1, pch = 16, bty = "n")
  
  #vsearch
  cols<-circleCol[factor(metaVSEARCH$ptInOut, levels = c("Outpatient", "Inpatient"))]
  
  MDS<-capscale(t(vsearch)~1,distance = "bray")
  percentVariance<-MDS$CA$eig/sum(eigenvals(MDS))*100
  pval<-adonis2(t(vsearch)~metaVSEARCH$ptInOut, method="bray")$`Pr(>F)`[1]
  statusPlot<-ordiplot(MDS,choices=c(1,2),type="none",cex.lab=1,
                       xlab=paste("MDS1  ", format(percentVariance[1], digits = 4), "%", sep = ""),
                       ylab=paste("MDS2  ", format(percentVariance[2], digits = 4), "%", sep = ""),
                       main=paste(timeFrames[i], "vsearch n =", ncol(vsearch), "\nP-value:",pval))
  points(statusPlot,"sites", pch=19, cex=2.5, col=adjustcolor(cols, alpha.f = 0.5))
  ordiellipse(statusPlot, metaVSEARCH$ptInOut, kind="se", conf=0.95, lwd=4, draw = "lines", col=circleCol[1],show.groups="Outpatient",label=T,font=2,cex=1) 
  ordiellipse(statusPlot, metaVSEARCH$ptInOut, kind="se", conf=0.95, lwd=4, draw = "lines", col=circleCol[2],show.groups="Inpatient",label=T,font=2,cex=1) 
  legend("topright", c("Outpatient", "Inpatient"), col = circleCol[1:2], cex = 1, pch = 16, bty = "n")
  
}

dev.off()


