#MDS plots for in and out patients.
rm(list = ls())
library("vegan")

metaData<-read.csv("metaWithBins.csv", header = T, row.names = 1)

#counts tables
brackenT<-read.csv("CountsTables/brackenFiltered.csv", header = T, row.names = 1, check.names = F)
amrT<-read.csv("CountsTables/amrFiltered.csv", header = T, row.names = 1, check.names = F)
rgiT<-read.csv("CountsTables/rgiFiltered.csv", header = T, row.names = 1, check.names = F)
vsearchT<-read.csv("CountsTables/vsearchFiltered.csv", header = T, row.names = 1, check.names = F)

metaBRACKEN<-metaData[colnames(brackenT), , drop = F]
metaAMR<-metaData[colnames(amrT), , drop = F]
metaRGI<-metaData[colnames(rgiT), , drop = F]
metaVSEARCH<-metaData[colnames(vsearchT), , drop = F]

pdf("Plots/MDSForInAndOutPatients.pdf", width=12, height=12)
par(mfrow=c(2,2))
par(mar=c(5,6,4,1)+.1)

#Bracken
circleCol<-c("cornflowerblue", "darkorange")
cols<-circleCol[factor(metaBRACKEN$ptInOut, levels = c("Outpatient", "Inpatient"))]

MDS<-capscale(t(brackenT)~1,distance = "bray")
percentVariance<-MDS$CA$eig/sum(eigenvals(MDS))*100
pval<-adonis2(t(brackenT)~metaBRACKEN$ptInOut, method="bray")$`Pr(>F)`[1]
statusPlot<-ordiplot(MDS,choices = c(1,2),type="none",cex.lab = 1,
                     xlab = paste("MDS1  ", format(percentVariance[1], digits = 4), "%", sep = ""),
                     ylab = paste("MDS2  ", format(percentVariance[2], digits = 4), "%", sep = ""),
                     main = paste("Bracken-Species n =", ncol(brackenT), "\nP-value:",pval))
points(statusPlot,"sites", pch = 19, cex = 2.5, col = adjustcolor(cols, alpha.f = 0.5))
ordiellipse(statusPlot, metaBRACKEN$ptInOut, kind="se", conf=0.95, lwd=4, draw = "lines", col=circleCol[1],show.groups="Outpatient",label=T,font=2,cex=1) 
ordiellipse(statusPlot, metaBRACKEN$ptInOut, kind="se", conf=0.95, lwd=4, draw = "lines", col=circleCol[2],show.groups="Inpatient",label=T,font=2,cex=1) 
legend("topright", c("Outpatient", "Inpatient"), col = circleCol[1:2], cex = 1, pch = 16, bty = "n")

#AMR
circleCol<-c("cornflowerblue", "darkorange")
cols<-circleCol[factor(metaAMR$ptInOut, levels = c("Outpatient", "Inpatient"))]

MDS<-capscale(t(amrT)~1,distance = "bray")
percentVariance<-MDS$CA$eig/sum(eigenvals(MDS))*100
pval<-adonis2(t(amrT)~metaAMR$ptInOut, method="bray")$`Pr(>F)`[1]
statusPlot<-ordiplot(MDS,choices=c(1,2),type="none",cex.lab=1,
                     xlab=paste("MDS1  ", format(percentVariance[1], digits = 4), "%", sep = ""),
                     ylab=paste("MDS2  ", format(percentVariance[2], digits = 4), "%", sep = ""),
                     main=paste("AMR n =", ncol(amrT), "\nP-value:",pval))
points(statusPlot,"sites", pch=19, cex=2.5, col=adjustcolor(cols, alpha.f = 0.5))
ordiellipse(statusPlot, metaAMR$ptInOut, kind="se", conf=0.95, lwd=4, draw = "lines", col=circleCol[1],show.groups="Outpatient",label=T,font=2,cex=1) 
ordiellipse(statusPlot, metaAMR$ptInOut, kind="se", conf=0.95, lwd=4, draw = "lines", col=circleCol[2],show.groups="Inpatient",label=T,font=2,cex=1) 
legend("topright", c("Outpatient", "Inpatient"), col = circleCol[1:2], cex = 1, pch = 16, bty = "n")

#RGI
circleCol<-c("cornflowerblue", "darkorange")
cols<-circleCol[factor(metaRGI$ptInOut, levels = c("Outpatient", "Inpatient"))]

MDS<-capscale(t(rgiT)~1,distance = "bray")
percentVariance<-MDS$CA$eig/sum(eigenvals(MDS))*100
pval<-adonis2(t(rgiT)~metaRGI$ptInOut, method="bray")$`Pr(>F)`[1]
statusPlot<-ordiplot(MDS,choices=c(1,2),type="none",cex.lab=1,
                     xlab=paste("MDS1  ", format(percentVariance[1], digits = 4), "%", sep = ""),
                     ylab=paste("MDS2  ", format(percentVariance[2], digits = 4), "%", sep = ""),
                     main=paste("RGI n =", ncol(rgiT), "\nP-value:",pval))
points(statusPlot,"sites", pch=19, cex=2.5, col=adjustcolor(cols, alpha.f = 0.5))
ordiellipse(statusPlot, metaRGI$ptInOut, kind="se", conf=0.95, lwd=4, draw = "lines", col=circleCol[1],show.groups="Outpatient",label=T,font=2,cex=1) 
ordiellipse(statusPlot, metaRGI$ptInOut, kind="se", conf=0.95, lwd=4, draw = "lines", col=circleCol[2],show.groups="Inpatient",label=T,font=2,cex=1) 
legend("topright", c("Outpatient", "Inpatient"), col = circleCol[1:2], cex = 1, pch = 16, bty = "n")

#vsearch
circleCol<-c("cornflowerblue", "darkorange")
cols<-circleCol[factor(metaVSEARCH$ptInOut, levels = c("Outpatient", "Inpatient"))]

MDS<-capscale(t(vsearchT)~1,distance = "bray")
percentVariance<-MDS$CA$eig/sum(eigenvals(MDS))*100
pval<-adonis2(t(vsearchT)~metaVSEARCH$ptInOut, method="bray")$`Pr(>F)`[1]
statusPlot<-ordiplot(MDS,choices=c(1,2),type="none",cex.lab=1,
                     xlab=paste("MDS1  ", format(percentVariance[1], digits = 4), "%", sep = ""),
                     ylab=paste("MDS2  ", format(percentVariance[2], digits = 4), "%", sep = ""),
                     main=paste("vsearch n =", ncol(vsearchT), "\nP-value:",pval))
points(statusPlot,"sites", pch=19, cex=2.5, col=adjustcolor(cols, alpha.f = 0.5))
ordiellipse(statusPlot, metaVSEARCH$ptInOut, kind="se", conf=0.95, lwd=4, draw = "lines", col=circleCol[1],show.groups="Outpatient",label=T,font=2,cex=1) 
ordiellipse(statusPlot, metaVSEARCH$ptInOut, kind="se", conf=0.95, lwd=4, draw = "lines", col=circleCol[2],show.groups="Inpatient",label=T,font=2,cex=1) 
legend("topright", c("Outpatient", "Inpatient"), col = circleCol[1:2], cex = 1, pch = 16, bty = "n")

dev.off()
