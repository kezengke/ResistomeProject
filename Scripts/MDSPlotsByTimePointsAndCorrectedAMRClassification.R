#MDS plots colored by timepoints based on MASTER_AMRlist_2023_03.csv Corrected AMR classification column
rm(list = ls())
library("vegan")
library(stringr)
library("dplyr")
library("RColorBrewer")

geneMatchT<-read.csv("MASTER_AMRlist_2023_03.csv", sep = ",", header = T)
metaData<-read.csv("metaWithBins.csv", header = T, row.names = 1)

#gene counts tables
amrT<-read.csv("CountsTables/amrNormalized.csv", header = T, row.names = 1, check.names = F)
rgiT<-read.csv("CountsTables/rgiNormalized.csv", header = T, row.names = 1, check.names = F)
vsearchT<-read.csv("CountsTables/vsearchNormalized.csv", header = T, row.names = 1, check.names = F)

#AMR
#keep only the rows present in amr counts table in master list
amrMatchT<-geneMatchT %>% filter(row_number() %in% match(rownames(amrT), geneMatchT$AMR.name))
amrMatchT<-na.omit(amrMatchT)
amrMatchT<-amrMatchT[, c(18,17)] #keep only the gene name and Corrected AMR classification column

amrTypes<-unique(amrMatchT$Corrected.AMR.classification)
newamrT<-c()
for (i in 1:length(amrTypes)) {
  genes<-amrMatchT %>% filter(amrMatchT$Corrected.AMR.classification == Corrected.AMR.classification[i])
  currentFrame<-amrT[genes$AMR.name, ]
  sumRow<-colSums(currentFrame)
  newamrT<-rbind(newamrT, sumRow)
}
rownames(newamrT)<-amrTypes
#exclude samples with 0 counts in every row
newamrT<-newamrT[, -which(colSums(newamrT) == 0)]

#RGI
#keep only the rows present in rgi counts table in master list
rgiMatchT<-geneMatchT %>% filter(row_number() %in% match(rownames(rgiT), geneMatchT$RGI.CARD.Short.Name))
rgiMatchT<-rgiMatchT[, c(15,17)] #keep only the gene name and Corrected AMR classification column

amrTypes<-unique(rgiMatchT$Corrected.AMR.classification)
newrgiT<-c()
for (i in 1:length(amrTypes)) {
  genes<-rgiMatchT %>% filter(rgiMatchT$Corrected.AMR.classification == amrTypes[i])
  currentFrame<-rgiT[genes$RGI.CARD.Short.Name, ]
  sumRow<-colSums(currentFrame)
  newrgiT<-rbind(newrgiT, sumRow)
}
rownames(newrgiT)<-amrTypes

#VSEARCH
#keep only the rows present in vsearch counts table in master list
vsearchMatchT<-geneMatchT %>% filter(row_number() %in% match(rownames(vsearchT), geneMatchT$Vsearch.ARO.Name))
vsearchMatchT<-vsearchMatchT[, c(7,17)] #keep only the gene name and Corrected AMR classification column

amrTypes<-vsearchMatchT$Corrected.AMR.classification
newvsearchT<-c()
for (i in 1:length(amrTypes)) {
  genes<-vsearchMatchT %>% filter(vsearchMatchT$Corrected.AMR.classification == amrTypes[i])
  currentFrame<-vsearchT[genes$Vsearch.ARO.Name, ]
  sumRow<-colSums(currentFrame)
  newvsearchT<-rbind(newvsearchT, sumRow)
}
rownames(newvsearchT)<-amrTypes

metaAMR<-metaData[colnames(newamrT), , drop = F]
metaRGI<-metaData[colnames(newrgiT), , drop = F]
metaVSEARCH<-metaData[colnames(newvsearchT), , drop = F]

pdf("Plots/MDSForCorrectedAMRClassification(ColoredByTimePoints).pdf", width=18, height=6)
par(mfrow=c(1,3))
par(mar=c(5,6,4,1)+.1)

circleCol<-brewer.pal(length(unique(metaAMR$bins)), "Spectral")
cols<-circleCol[factor(metaAMR$bins, levels = c("PRE", "D0", "D7", "D14", "D21", "D28", "D35", "D60", "D100", "D180"))]
MDS<-capscale(t(newamrT)~1,distance = "bray")
percentVariance<-MDS$CA$eig/sum(eigenvals(MDS))*100
pval<-adonis2(t(newamrT)~metaAMR$bins, method="bray")$`Pr(>F)`[1]
statusPlot<-ordiplot(MDS,choices=c(1,2),type="none",cex.lab=1,
                     xlab=paste("MDS1  ", format(percentVariance[1], digits = 4), "%", sep = ""),
                     ylab=paste("MDS2  ", format(percentVariance[2], digits = 4), "%", sep = ""),
                     main=paste("AMR(CorrectedAMRClassification) n =", ncol(newamrT), "\nP-value:",pval))
points(statusPlot,"sites", pch=19, cex=2.5, col=adjustcolor(cols, alpha.f = 0.5))
ordiellipse(statusPlot, metaAMR$bins, kind="se", conf=0.95, lwd=4, draw = "lines", col=circleCol) 
legend("topleft", c("PRE", "D0", "D7", "D14", "D21", "D28", "D35", "D60", "D100", "D180"), col = circleCol[1:10], cex = 1, pch = 16, bty = "n")

circleCol<-brewer.pal(length(unique(metaRGI$bins)), "Spectral")
cols<-circleCol[factor(metaRGI$bins, levels = c("PRE", "D0", "D7", "D14", "D21", "D28", "D35", "D60", "D100", "D180"))]
MDS<-capscale(t(newrgiT)~1,distance = "bray")
percentVariance<-MDS$CA$eig/sum(eigenvals(MDS))*100
pval<-adonis2(t(newrgiT)~metaRGI$bins, method="bray")$`Pr(>F)`[1]
statusPlot<-ordiplot(MDS,choices=c(1,2),type="none",cex.lab=1,
                     xlab=paste("MDS1  ", format(percentVariance[1], digits = 4), "%", sep = ""),
                     ylab=paste("MDS2  ", format(percentVariance[2], digits = 4), "%", sep = ""),
                     main=paste("RGI(CorrectedAMRClassification) n =", ncol(newrgiT), "\nP-value:",pval))
points(statusPlot,"sites", pch=19, cex=2.5, col=adjustcolor(cols, alpha.f = 0.5))
ordiellipse(statusPlot, metaRGI$bins, kind="se", conf=0.95, lwd=4, draw = "lines", col=circleCol) 
legend("topleft", c("PRE", "D0", "D7", "D14", "D21", "D28", "D35", "D60", "D100", "D180"), col = circleCol[1:10], cex = 1, pch = 16, bty = "n")

circleCol<-brewer.pal(length(unique(metaVSEARCH$bins)), "Spectral")
cols<-circleCol[factor(metaVSEARCH$bins, levels = c("PRE", "D0", "D7", "D14", "D21", "D28", "D35", "D60", "D100", "D180"))]
MDS<-capscale(t(newvsearchT)~1,distance = "bray")
percentVariance<-MDS$CA$eig/sum(eigenvals(MDS))*100
pval<-adonis2(t(newvsearchT)~metaVSEARCH$bins, method="bray")$`Pr(>F)`[1]
statusPlot<-ordiplot(MDS,choices=c(1,2),type="none",cex.lab=1,
                     xlab=paste("MDS1  ", format(percentVariance[1], digits = 4), "%", sep = ""),
                     ylab=paste("MDS2  ", format(percentVariance[2], digits = 4), "%", sep = ""),
                     main=paste("vsearch(CorrectedAMRClassification) n =", ncol(newvsearchT), "\nP-value:",pval))
points(statusPlot,"sites", pch=19, cex=2.5, col=adjustcolor(cols, alpha.f = 0.5))
ordiellipse(statusPlot, metaVSEARCH$bins, kind="se", conf=0.95, lwd=4, draw = "lines", col=circleCol) 
legend("topright", c("PRE", "D0", "D7", "D14", "D21", "D28", "D35", "D60", "D100", "D180"), col = circleCol[1:10], cex = 1, pch = 16, bty = "n")

dev.off()
