#LM for MDS1 and MDS2
rm(list = ls())
library("vegan")
library("nlme")
#metadata

metaData<-read.csv("metaWithBins.csv", header = T, row.names = 1)

#gene counts tables
speciesT<-read.csv("CountsTables/brackenFiltered.csv", header = T, row.names = 1, check.names = F)
genusT<-read.csv("CountsTables/genusFiltered.csv", header = T, row.names = 1, check.names = F)
amrT<-read.csv("CountsTables/amrFiltered.csv", header = T, row.names = 1, check.names = F)
rgiT<-read.csv("CountsTables/rgiFiltered.csv", header = T, row.names = 1, check.names = F)
vsearchT<-read.csv("CountsTables/vsearchFiltered.csv", header = T, row.names = 1, check.names = F)

metaSPECIES<-metaData[colnames(speciesT), , drop = F]
metaGENUS<-metaData[colnames(genusT), , drop = F]
metaAMR<-metaData[colnames(amrT), , drop = F]
metaRGI<-metaData[colnames(rgiT), , drop = F]
metaVSEARCH<-metaData[colnames(vsearchT), , drop = F]

plotFun <- function(mdsNum, MDS, meta, Tname) {
  myM<-data.frame(MDS[, mdsNum], meta$Timepoint, meta$ptInOut, meta$ID)
  colnames(myM)<-c("counts", "timePoint", "InOut", "ID")
  Model<-lme(counts ~ timePoint * InOut, random = ~1 | ID, data = myM)
  pvals<-summary(Model)$tTable[c(2:4), 5]
  group_colors<-c("cornflowerblue", "darkorange")[as.numeric(as.factor(meta$ptInOut))]
  plot(log10(meta$Timepoint), MDS[, mdsNum], 
       pch = 19, col = group_colors, 
       main = paste0("MDS", mdsNum, "(", Tname, ")", 
                    "\nTimePoint(P):", signif(pvals[1], 3),
                    "\nInOut(P):", signif(pvals[2], 3),
                    "\nInteraction(P):", signif( pvals[2], 3)), 
       cex.main = 1, xlab = "log10(TimePoints)", ylab = paste0("MDS", mdsNum))
  lines(log10(sort(meta$Timepoint)), fitted(Model)[order(meta$Timepoint)], 
        col = "grey", type = "l")
  legend("bottomright", c("Outpatient", "Inpatient"), col = c("cornflowerblue", "darkorange"), cex = 1, pch = 16, bty = "n")
}

pdf("Plots/MDSlms.pdf", width=12, height=18)
par(mfrow=c(3, 2))
par(mar=c(5, 6, 4, 1)+.1)
#Species
MDS<-summary(capscale(t(speciesT)~1,distance = "bray"))[["sites"]][,c(1,2,3,4,5)]

plotFun(1, MDS, metaSPECIES, "Species")
plotFun(2, MDS, metaSPECIES, "Species")

#Genus
MDS<-summary(capscale(t(genusT)~1,distance = "bray"))[["sites"]][,c(1,2,3,4,5)]

plotFun(1, MDS, metaGENUS, "Genus")
plotFun(2, MDS, metaGENUS, "Genus")

#AMR
MDS<-summary(capscale(t(amrT)~1,distance = "bray"))[["sites"]][,c(1,2,3,4,5)]

plotFun(1, MDS, metaAMR, "AMR")
plotFun(2, MDS, metaAMR, "AMR")

#RGI
MDS<-summary(capscale(t(rgiT)~1,distance = "bray"))[["sites"]][,c(1,2,3,4,5)]

plotFun(1, MDS, metaRGI, "RGI")
plotFun(2, MDS, metaRGI, "RGI")

#vsearch
MDS<-summary(capscale(t(vsearchT)~1,distance = "bray"))[["sites"]][,c(1,2,3,4,5)]

plotFun(1, MDS, metaVSEARCH, "VSEARCH")
plotFun(2, MDS, metaVSEARCH, "VSEARCH")

dev.off()

