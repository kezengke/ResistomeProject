rm(list = ls())
library(stringr)
library("RColorBrewer")
library("vegan")
library("dplyr")

metaData<-read.csv("metaWithBins.csv", header = T, row.names = 1)

#counts tables
brackenT<-read.csv("CountsTables/brackenNormalized.csv", header = T, row.names = 1, check.names = F)
amrT<-read.csv("CountsTables/amrNormalized.csv", header = T, row.names = 1, check.names = F)
rgiT<-read.csv("CountsTables/rgiNormalized.csv", header = T, row.names = 1, check.names = F)
vsearchT<-read.csv("CountsTables/vsearchNormalized.csv", header = T, row.names = 1, check.names = F)

#exclude days after 600 days
brackenT<-brackenT[, -which(colnames(brackenT) %in% c(rownames(metaData)[metaData$Timepoint>600])) ]
amrT<-amrT[, -which(colnames(amrT) %in% c(rownames(metaData)[metaData$Timepoint>600])) ]
rgiT<-rgiT[, -which(colnames(rgiT) %in% c(rownames(metaData)[metaData$Timepoint>600])) ]
vsearchT<-vsearchT[, -which(colnames(vsearchT) %in% c(rownames(metaData)[metaData$Timepoint>600])) ]

#meta for each table
metaBRACKEN<-metaData[colnames(brackenT), , drop = F]
metaAMR<-metaData[colnames(amrT), , drop = F]
metaRGI<-metaData[colnames(rgiT), , drop = F]
metaVSEARCH<-metaData[colnames(vsearchT), , drop = F]

pdf("Plots/MDSForAllLessThan600Days(ColoredByTimePoints).pdf", width=12, height=12)
par(mfrow=c(2,2))
par(mar=c(5,6,4,1)+.1)

#Bracken
circleCol<-brewer.pal(length(unique(metaBRACKEN$bins)), "Spectral")
cols<-circleCol[factor(metaBRACKEN$bins, levels = c("PRE", "D0", "D7", "D14", "D21", "D28", "D35", "D60", "D100", "D180"))]

MDS<-capscale(t(brackenT)~1,distance = "bray")
percentVariance<-MDS$CA$eig/sum(eigenvals(MDS))*100
pval<-adonis2(t(brackenT)~metaBRACKEN$bins, method="bray")$`Pr(>F)`[1]
statusPlot<-ordiplot(MDS,choices = c(1,2),type="none",cex.lab = 1,
                     xlab = paste("MDS1  ", format(percentVariance[1], digits = 4), "%", sep = ""),
                     ylab = paste("MDS2  ", format(percentVariance[2], digits = 4), "%", sep = ""),
                     main = paste("Bracken-Species n =", ncol(brackenT), "\nP-value:",pval))
points(statusPlot,"sites", pch = 19, cex = 2.5, col = adjustcolor(cols, alpha.f = 0.5))
for (i in 1:length(unique(metaBRACKEN$bins))) {
  ordiellipse(statusPlot, metaBRACKEN$bins, kind="se", conf=0.95, lwd=4, draw = "lines", col=circleCol[i], 
              show.groups = c("PRE", "D0", "D7", "D14", "D21", "D28", "D35", "D60", "D100", "D180")[i],
              font=2,cex=1) 
}
legend("topright", c("PRE", "D0", "D7", "D14", "D21", "D28", "D35", "D60", "D100", "D180"), col = circleCol[1:10], cex = 1, pch = 16, bty = "n")

#AMR
circleCol<-brewer.pal(length(unique(metaAMR$bins)), "Spectral")
cols<-circleCol[factor(metaAMR$bins, levels = c("PRE", "D0", "D7", "D14", "D21", "D28", "D35", "D60", "D100", "D180"))]

MDS<-capscale(t(amrT)~1,distance = "bray")
percentVariance<-MDS$CA$eig/sum(eigenvals(MDS))*100
pval<-adonis2(t(amrT)~metaAMR$bins, method="bray")$`Pr(>F)`[1]
statusPlot<-ordiplot(MDS,choices=c(1,2),type="none",cex.lab=1,
                     xlab=paste("MDS1  ", format(percentVariance[1], digits = 4), "%", sep = ""),
                     ylab=paste("MDS2  ", format(percentVariance[2], digits = 4), "%", sep = ""),
                     main=paste("AMR n =", ncol(amrT), "\nP-value:",pval))
points(statusPlot,"sites", pch=19, cex=2.5, col=adjustcolor(cols, alpha.f = 0.5))
for (i in 1:length(unique(metaAMR$bins))) {
  ordiellipse(statusPlot, metaAMR$bins, kind="se", conf=0.95, lwd=4, draw = "lines", col=circleCol[i], 
              show.groups = c("PRE", "D0", "D7", "D14", "D21", "D28", "D35", "D60", "D100", "D180")[i],
              font=2,cex=1) 
}
legend("topleft", c("PRE", "D0", "D7", "D14", "D21", "D28", "D35", "D60", "D100", "D180"), col = circleCol[1:10], cex = 1, pch = 16, bty = "n")

#RGI
circleCol<-brewer.pal(length(unique(metaRGI$bins)), "Spectral")
cols<-circleCol[factor(metaRGI$bins, levels = c("PRE", "D0", "D7", "D14", "D21", "D28", "D35", "D60", "D100", "D180"))]

MDS<-capscale(t(rgiT)~1,distance = "bray")
percentVariance<-MDS$CA$eig/sum(eigenvals(MDS))*100
pval<-adonis2(t(rgiT)~metaRGI$bins, method="bray")$`Pr(>F)`[1]
statusPlot<-ordiplot(MDS,choices=c(1,2),type="none",cex.lab=1,
                     xlab=paste("MDS1  ", format(percentVariance[1], digits = 4), "%", sep = ""),
                     ylab=paste("MDS2  ", format(percentVariance[2], digits = 4), "%", sep = ""),
                     main=paste("RGI n =", ncol(rgiT), "\nP-value:",pval))
points(statusPlot,"sites", pch=19, cex=2.5, col=adjustcolor(cols, alpha.f = 0.5))
for (i in 1:length(unique(metaRGI$bins))) {
  ordiellipse(statusPlot, metaRGI$bins, kind="se", conf=0.95, lwd=4, draw = "lines", col=circleCol[i], 
              show.groups = c("PRE", "D0", "D7", "D14", "D21", "D28", "D35", "D60", "D100", "D180")[i],
              font=2,cex=1) 
}
legend("topright", c("PRE", "D0", "D7", "D14", "D21", "D28", "D35", "D60", "D100", "D180"), col = circleCol[1:10], cex = 1, pch = 16, bty = "n")

#vsearch
circleCol<-brewer.pal(length(unique(metaVSEARCH$bins)), "Spectral") 
cols<-circleCol[factor(metaVSEARCH$bins, levels = c("PRE", "D0", "D7", "D14", "D21", "D28", "D35", "D60", "D100", "D180"))]

MDS<-capscale(t(vsearchT)~1,distance = "bray")
percentVariance<-MDS$CA$eig/sum(eigenvals(MDS))*100
pval<-adonis2(t(vsearchT)~metaVSEARCH$bins, method="bray")$`Pr(>F)`[1]
statusPlot<-ordiplot(MDS,choices=c(1,2),type="none",cex.lab=1,
                     xlab=paste("MDS1  ", format(percentVariance[1], digits = 4), "%", sep = ""),
                     ylab=paste("MDS2  ", format(percentVariance[2], digits = 4), "%", sep = ""),
                     main=paste("vsearch n =", ncol(vsearchT), "\nP-value:",pval))
points(statusPlot,"sites", pch=19, cex=2.5, col=adjustcolor(cols, alpha.f = 0.5))
for (i in 1:length(unique(metaVSEARCH$bins))) {
  ordiellipse(statusPlot, metaVSEARCH$bins, kind="se", conf=0.95, lwd=4, draw = "lines", col=circleCol[i], 
              show.groups = c("PRE", "D0", "D7", "D14", "D21", "D28", "D35", "D60", "D100", "D180")[i],
              font=2,cex=1) 
}
legend("topright", c("PRE", "D0", "D7", "D14", "D21", "D28", "D35", "D60", "D100", "D180"), col = circleCol[1:10], cex = 1, pch = 16, bty = "n")

dev.off()
