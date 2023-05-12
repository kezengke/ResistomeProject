#MDS plots for in and out patients for pre and post.
rm(list = ls())
library("vegan")

metaData<-read.csv("metaWithBins.csv", header = T, row.names = 1)

#counts tables
brackenT<-read.csv("CountsTables/brackenNormalized.csv", header = T, row.names = 1, check.names = F)
amrT<-read.csv("CountsTables/amrNormalized.csv", header = T, row.names = 1, check.names = F)
rgiT<-read.csv("CountsTables/rgiNormalized.csv", header = T, row.names = 1, check.names = F)
vsearchT<-read.csv("CountsTables/vsearchNormalized.csv", header = T, row.names = 1, check.names = F)

pres<-rownames(metaData)[metaData$bins == "PRE"]
posts<-rownames(metaData)[metaData$bins != "PRE"]

brackenPre<-brackenT[, pres, drop = F]
brackenPost<-brackenT[, posts, drop = F]
amrPre<-amrT[, pres, drop = F]
amrPost<-amrT[, intersect(colnames(amrT), posts), drop = F]
rgiPre<-rgiT[, pres, drop = F]
rgiPost<-rgiT[, posts, drop = F]
vsearchPre<-vsearchT[, pres, drop = F]
vsearchPost<-vsearchT[, posts, drop = F]

metaBRACKENpre<-metaData[colnames(brackenPre), , drop = F]
metaBRACKENpost<-metaData[colnames(brackenPost), , drop = F]
metaAMRpre<-metaData[colnames(amrPre), , drop = F]
metaAMRpost<-metaData[colnames(amrPost), , drop = F]
metaRGIpre<-metaData[colnames(rgiPre), , drop = F]
metaRGIpost<-metaData[colnames(rgiPost), , drop = F]
metaVSEARCHpre<-metaData[colnames(vsearchPre), , drop = F]
metaVSEARCHpost<-metaData[colnames(vsearchPost), , drop = F]

pdf("Plots/MDSForInAndOutPatients(PrePost).pdf", width=24, height=12)
par(mfrow=c(2,4))
par(mar=c(5,6,4,1)+.1)

#####PRE#####
#Bracken
circleCol<-c("cornflowerblue", "darkorange")
cols<-circleCol[factor(metaBRACKENpre$ptInOut, levels = c("Outpatient", "Inpatient"))]

MDS<-capscale(t(brackenPre)~1,distance = "bray")
percentVariance<-MDS$CA$eig/sum(eigenvals(MDS))*100
pval<-adonis2(t(brackenPre)~metaBRACKENpre$ptInOut, method="bray")$`Pr(>F)`[1]
statusPlot<-ordiplot(MDS,choices = c(1,2),type="none",cex.lab = 1,
                     xlab = paste("MDS1  ", format(percentVariance[1], digits = 4), "%", sep = ""),
                     ylab = paste("MDS2  ", format(percentVariance[2], digits = 4), "%", sep = ""),
                     main = paste("(PRE)Bracken-Species n =", ncol(brackenPre), "\nP-value:",pval))
points(statusPlot,"sites", pch = 19, cex = 2.5, col = adjustcolor(cols, alpha.f = 0.5))
ordiellipse(statusPlot, metaBRACKENpre$ptInOut, kind="se", conf=0.95, lwd=4, draw = "lines", col=circleCol[1],show.groups="Outpatient",label=T,font=2,cex=1) 
ordiellipse(statusPlot, metaBRACKENpre$ptInOut, kind="se", conf=0.95, lwd=4, draw = "lines", col=circleCol[2],show.groups="Inpatient",label=T,font=2,cex=1) 
legend("topright", c("Outpatient", "Inpatient"), col = circleCol[1:2], cex = 1, pch = 16, bty = "n")

#AMR
circleCol<-c("cornflowerblue", "darkorange")
cols<-circleCol[factor(metaAMRpre$ptInOut, levels = c("Outpatient", "Inpatient"))]

MDS<-capscale(t(amrPre)~1,distance = "bray")
percentVariance<-MDS$CA$eig/sum(eigenvals(MDS))*100
pval<-adonis2(t(amrPre)~metaAMRpre$ptInOut, method="bray")$`Pr(>F)`[1]
statusPlot<-ordiplot(MDS,choices=c(1,2),type="none",cex.lab=1,
                     xlab=paste("MDS1  ", format(percentVariance[1], digits = 4), "%", sep = ""),
                     ylab=paste("MDS2  ", format(percentVariance[2], digits = 4), "%", sep = ""),
                     main=paste("(PRE)AMR n =", ncol(amrPre), "\nP-value:",pval))
points(statusPlot,"sites", pch=19, cex=2.5, col=adjustcolor(cols, alpha.f = 0.5))
ordiellipse(statusPlot, metaAMRpre$ptInOut, kind="se", conf=0.95, lwd=4, draw = "lines", col=circleCol[1],show.groups="Outpatient",label=T,font=2,cex=1) 
ordiellipse(statusPlot, metaAMRpre$ptInOut, kind="se", conf=0.95, lwd=4, draw = "lines", col=circleCol[2],show.groups="Inpatient",label=T,font=2,cex=1) 
legend("topright", c("Outpatient", "Inpatient"), col = circleCol[1:2], cex = 1, pch = 16, bty = "n")

#RGI
circleCol<-c("cornflowerblue", "darkorange")
cols<-circleCol[factor(metaRGIpre$ptInOut, levels = c("Outpatient", "Inpatient"))]

MDS<-capscale(t(rgiPre)~1,distance = "bray")
percentVariance<-MDS$CA$eig/sum(eigenvals(MDS))*100
pval<-adonis2(t(rgiPre)~metaRGIpre$ptInOut, method="bray")$`Pr(>F)`[1]
statusPlot<-ordiplot(MDS,choices=c(1,2),type="none",cex.lab=1,
                     xlab=paste("MDS1  ", format(percentVariance[1], digits = 4), "%", sep = ""),
                     ylab=paste("MDS2  ", format(percentVariance[2], digits = 4), "%", sep = ""),
                     main=paste("(PRE)RGI n =", ncol(rgiPre), "\nP-value:",pval))
points(statusPlot,"sites", pch=19, cex=2.5, col=adjustcolor(cols, alpha.f = 0.5))
ordiellipse(statusPlot, metaRGIpre$ptInOut, kind="se", conf=0.95, lwd=4, draw = "lines", col=circleCol[1],show.groups="Outpatient",label=T,font=2,cex=1) 
ordiellipse(statusPlot, metaRGIpre$ptInOut, kind="se", conf=0.95, lwd=4, draw = "lines", col=circleCol[2],show.groups="Inpatient",label=T,font=2,cex=1) 
legend("topright", c("Outpatient", "Inpatient"), col = circleCol[1:2], cex = 1, pch = 16, bty = "n")

#vsearch
circleCol<-c("cornflowerblue", "darkorange")
cols<-circleCol[factor(metaVSEARCHpre$ptInOut, levels = c("Outpatient", "Inpatient"))]

MDS<-capscale(t(vsearchPre)~1,distance = "bray")
percentVariance<-MDS$CA$eig/sum(eigenvals(MDS))*100
pval<-adonis2(t(vsearchPre)~metaVSEARCHpre$ptInOut, method="bray")$`Pr(>F)`[1]
statusPlot<-ordiplot(MDS,choices=c(1,2),type="none",cex.lab=1,
                     xlab=paste("MDS1  ", format(percentVariance[1], digits = 4), "%", sep = ""),
                     ylab=paste("MDS2  ", format(percentVariance[2], digits = 4), "%", sep = ""),
                     main=paste("(PRE)vsearch n =", ncol(vsearchPre), "\nP-value:",pval))
points(statusPlot,"sites", pch=19, cex=2.5, col=adjustcolor(cols, alpha.f = 0.5))
ordiellipse(statusPlot, metaVSEARCHpre$ptInOut, kind="se", conf=0.95, lwd=4, draw = "lines", col=circleCol[1],show.groups="Outpatient",label=T,font=2,cex=1) 
ordiellipse(statusPlot, metaVSEARCHpre$ptInOut, kind="se", conf=0.95, lwd=4, draw = "lines", col=circleCol[2],show.groups="Inpatient",label=T,font=2,cex=1) 
legend("topright", c("Outpatient", "Inpatient"), col = circleCol[1:2], cex = 1, pch = 16, bty = "n")

#####POST#####
#Bracken
circleCol<-c("cornflowerblue", "darkorange")
cols<-circleCol[factor(metaBRACKENpost$ptInOut, levels = c("Outpatient", "Inpatient"))]

MDS<-capscale(t(brackenPost)~1,distance = "bray")
percentVariance<-MDS$CA$eig/sum(eigenvals(MDS))*100
pval<-adonis2(t(brackenPost)~metaBRACKENpost$ptInOut, method="bray")$`Pr(>F)`[1]
statusPlot<-ordiplot(MDS,choices = c(1,2),type="none",cex.lab = 1,
                     xlab = paste("MDS1  ", format(percentVariance[1], digits = 4), "%", sep = ""),
                     ylab = paste("MDS2  ", format(percentVariance[2], digits = 4), "%", sep = ""),
                     main = paste("(POST)Bracken-Species n =", ncol(brackenPost), "\nP-value:",pval))
points(statusPlot,"sites", pch = 19, cex = 2.5, col = adjustcolor(cols, alpha.f = 0.5))
ordiellipse(statusPlot, metaBRACKENpost$ptInOut, kind="se", conf=0.95, lwd=4, draw = "lines", col=circleCol[1],show.groups="Outpatient",label=T,font=2,cex=1) 
ordiellipse(statusPlot, metaBRACKENpost$ptInOut, kind="se", conf=0.95, lwd=4, draw = "lines", col=circleCol[2],show.groups="Inpatient",label=T,font=2,cex=1) 
legend("topright", c("Outpatient", "Inpatient"), col = circleCol[1:2], cex = 1, pch = 16, bty = "n")

#AMR
circleCol<-c("cornflowerblue", "darkorange")
cols<-circleCol[factor(metaAMRpost$ptInOut, levels = c("Outpatient", "Inpatient"))]

MDS<-capscale(t(amrPost)~1,distance = "bray")
percentVariance<-MDS$CA$eig/sum(eigenvals(MDS))*100
pval<-adonis2(t(amrPost)~metaAMRpost$ptInOut, method="bray")$`Pr(>F)`[1]
statusPlot<-ordiplot(MDS,choices=c(1,2),type="none",cex.lab=1,
                     xlab=paste("MDS1  ", format(percentVariance[1], digits = 4), "%", sep = ""),
                     ylab=paste("MDS2  ", format(percentVariance[2], digits = 4), "%", sep = ""),
                     main=paste("(POST)AMR n =", ncol(amrPost), "\nP-value:",pval))
points(statusPlot,"sites", pch=19, cex=2.5, col=adjustcolor(cols, alpha.f = 0.5))
ordiellipse(statusPlot, metaAMRpost$ptInOut, kind="se", conf=0.95, lwd=4, draw = "lines", col=circleCol[1],show.groups="Outpatient",label=T,font=2,cex=1) 
ordiellipse(statusPlot, metaAMRpost$ptInOut, kind="se", conf=0.95, lwd=4, draw = "lines", col=circleCol[2],show.groups="Inpatient",label=T,font=2,cex=1) 
legend("topright", c("Outpatient", "Inpatient"), col = circleCol[1:2], cex = 1, pch = 16, bty = "n")

#RGI
circleCol<-c("cornflowerblue", "darkorange")
cols<-circleCol[factor(metaRGIpost$ptInOut, levels = c("Outpatient", "Inpatient"))]

MDS<-capscale(t(rgiPost)~1,distance = "bray")
percentVariance<-MDS$CA$eig/sum(eigenvals(MDS))*100
pval<-adonis2(t(rgiPost)~metaRGIpost$ptInOut, method="bray")$`Pr(>F)`[1]
statusPlot<-ordiplot(MDS,choices=c(1,2),type="none",cex.lab=1,
                     xlab=paste("MDS1  ", format(percentVariance[1], digits = 4), "%", sep = ""),
                     ylab=paste("MDS2  ", format(percentVariance[2], digits = 4), "%", sep = ""),
                     main=paste("(POST)RGI n =", ncol(rgiPost), "\nP-value:",pval))
points(statusPlot,"sites", pch=19, cex=2.5, col=adjustcolor(cols, alpha.f = 0.5))
ordiellipse(statusPlot, metaRGIpost$ptInOut, kind="se", conf=0.95, lwd=4, draw = "lines", col=circleCol[1],show.groups="Outpatient",label=T,font=2,cex=1) 
ordiellipse(statusPlot, metaRGIpost$ptInOut, kind="se", conf=0.95, lwd=4, draw = "lines", col=circleCol[2],show.groups="Inpatient",label=T,font=2,cex=1) 
legend("topright", c("Outpatient", "Inpatient"), col = circleCol[1:2], cex = 1, pch = 16, bty = "n")

#vsearch
circleCol<-c("cornflowerblue", "darkorange")
cols<-circleCol[factor(metaVSEARCHpost$ptInOut, levels = c("Outpatient", "Inpatient"))]

MDS<-capscale(t(vsearchPost)~1,distance = "bray")
percentVariance<-MDS$CA$eig/sum(eigenvals(MDS))*100
pval<-adonis2(t(vsearchPost)~metaVSEARCHpost$ptInOut, method="bray")$`Pr(>F)`[1]
statusPlot<-ordiplot(MDS,choices=c(1,2),type="none",cex.lab=1,
                     xlab=paste("MDS1  ", format(percentVariance[1], digits = 4), "%", sep = ""),
                     ylab=paste("MDS2  ", format(percentVariance[2], digits = 4), "%", sep = ""),
                     main=paste("(POST)vsearch n =", ncol(vsearchPost), "\nP-value:",pval))
points(statusPlot,"sites", pch=19, cex=2.5, col=adjustcolor(cols, alpha.f = 0.5))
ordiellipse(statusPlot, metaVSEARCHpost$ptInOut, kind="se", conf=0.95, lwd=4, draw = "lines", col=circleCol[1],show.groups="Outpatient",label=T,font=2,cex=1) 
ordiellipse(statusPlot, metaVSEARCHpost$ptInOut, kind="se", conf=0.95, lwd=4, draw = "lines", col=circleCol[2],show.groups="Inpatient",label=T,font=2,cex=1) 
legend("topright", c("Outpatient", "Inpatient"), col = circleCol[1:2], cex = 1, pch = 16, bty = "n")

dev.off()
