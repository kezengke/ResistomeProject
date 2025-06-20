rm(list = ls())
library(stringr)
library("RColorBrewer")
library("vegan")

#metadata
metaData<-read.csv("metaWithBins.csv", header = T, row.names = 1)

#gene counts tables
brackenT<-read.csv("CountsTables/brackenFiltered.csv", header = T, row.names = 1, check.names = F)
amrT<-read.csv("CountsTables/amrFiltered.csv", header = T, row.names = 1, check.names = F)
rgiT<-read.csv("CountsTables/rgiFiltered.csv", header = T, row.names = 1, check.names = F)
vsearchT<-read.csv("CountsTables/vsearchFiltered.csv", header = T, row.names = 1, check.names = F)

metaBRACKEN<-metaData[colnames(brackenT), , drop = F]
metaAMR<-metaData[colnames(amrT), , drop = F]
metaRGI<-metaData[colnames(rgiT), , drop = F]
metaVSEARCH<-metaData[colnames(vsearchT), , drop = F]

pdf("Plots/MDS1And2BoxPlot(ByPatients).pdf", width=12, height=24)
par(mfrow=c(4,2))
par(mar=c(5,6,4,1)+.1)

#Bracken
MDS<-summary(capscale(t(brackenT)~1,distance = "bray"))[["sites"]][,c(1,2)]
pval1<-summary(aov(MDS[,1]~metaBRACKEN$ID))[[1]][["Pr(>F)"]][1]
pval2<-summary(aov(MDS[,2]~metaBRACKEN$ID))[[1]][["Pr(>F)"]][1]
boxplot(MDS[,1]~metaBRACKEN$ID, outline = F, col = "tan2", las = 3,
        ylab = "MDS1", xlab = "", main = paste0("Bracken-Species \nANOVA pvalue:", signif(pval1, digits = 3)))
boxplot(MDS[,2]~metaBRACKEN$ID, outline = F, col = "tan2", las = 3,
        ylab = "MDS2", xlab = "", main = paste0("Bracken-Species \nANOVA pvalue:", signif(pval2, digits = 3)))

#AMR
MDS<-summary(capscale(t(amrT)~1,distance = "bray"))[["sites"]][,c(1,2)]
pval1<-summary(aov(MDS[,1]~metaAMR$ID))[[1]][["Pr(>F)"]][1]
pval2<-summary(aov(MDS[,2]~metaAMR$ID))[[1]][["Pr(>F)"]][1]
boxplot(MDS[,1]~metaAMR$ID, outline = F, col = "coral3", las = 3,
        ylab = "MDS1", xlab = "", main = paste0("AMR \nANOVA pvalue:", signif(pval1, digits = 3)))
boxplot(MDS[,2]~metaAMR$ID, outline = F, col = "coral3", las = 3,
        ylab = "MDS2", xlab = "", main = paste0("AMR \nANOVA pvalue:", signif(pval2, digits = 3)))

#RGI
MDS<-summary(capscale(t(rgiT)~1,distance = "bray"))[["sites"]][,c(1,2)]
pval1<-summary(aov(MDS[,1]~metaRGI$ID))[[1]][["Pr(>F)"]][1]
pval2<-summary(aov(MDS[,2]~metaRGI$ID))[[1]][["Pr(>F)"]][1]
boxplot(MDS[,1]~metaRGI$ID, outline = F, col = "cornflowerblue", las = 3,
        ylab = "MDS1", xlab = "", main = paste0("RGI \nANOVA pvalue:", signif(pval1, digits = 3)))
boxplot(MDS[,2]~metaRGI$ID, outline = F, col = "cornflowerblue", las = 3,
        ylab = "MDS2", xlab = "", main = paste0("RGI \nANOVA pvalue:", signif(pval2, digits = 3)))

#vsearch
MDS<-summary(capscale(t(vsearchT)~1,distance = "bray"))[["sites"]][,c(1,2)]
pval1<-summary(aov(MDS[,1]~metaVSEARCH$ID))[[1]][["Pr(>F)"]][1]
pval2<-summary(aov(MDS[,1]~metaVSEARCH$ID))[[1]][["Pr(>F)"]][1]
boxplot(MDS[,1]~metaVSEARCH$ID, outline = F, col = "olivedrab4", las = 3,
        ylab = "MDS1", xlab = "", main = paste0("vsearch \nANOVA pvalue:",signif(pval1, digits = 3)))
boxplot(MDS[,2]~metaVSEARCH$ID, outline = F, col = "olivedrab4", las = 3,
        ylab = "MDS2", xlab = "", main = paste0("vsearch \nANOVA pvalue:",signif(pval2, digits = 3)))

dev.off()
