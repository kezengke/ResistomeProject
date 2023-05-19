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

#Bracken
pdf("Plots/TimeVsMDS1234(Bracken).pdf", width=12, height=12)
par(mfrow=c(2,2))
par(mar=c(5,6,4,1)+.1)
MDS<-summary(capscale(t(brackenT)~1,distance = "bray"))[["sites"]][,c(1,2,3,4)]
MDS<-log10(MDS+10)
plot(metaBRACKEN[rownames(MDS), 2], MDS[ ,1], 
     col = "tan2", pch = 19, 
     xlab = "Days", ylab = "log10(MDS1+10)",
     main = "Bracken MDS1")

plot(metaBRACKEN[rownames(MDS), 2], MDS[ ,2], 
     col = "tan2", pch = 19, 
     xlab = "Days", ylab = "log10(MDS2+10)",
     main = "Bracken MDS2")

plot(metaBRACKEN[rownames(MDS), 2], MDS[ ,3], 
     col = "tan2", pch = 19, 
     xlab = "Days", ylab = "log10(MDS3+10)",
     main = "Bracken MDS3")

plot(metaBRACKEN[rownames(MDS), 2], MDS[ ,4], 
     col = "tan2", pch = 19, 
     xlab = "Days", ylab = "log10(MDS4+10)",
     main = "Bracken MDS4")
dev.off()

#AMR
pdf("Plots/TimeVsMDS1234(AMR).pdf", width=12, height=12)
par(mfrow=c(2,2))
par(mar=c(5,6,4,1)+.1)
MDS<-summary(capscale(t(amrT)~1,distance = "bray"))[["sites"]][,c(1,2,3,4)]
MDS<-log10(MDS+10)
plot(metaAMR[rownames(MDS), 2], MDS[ ,1], 
     col = "coral3", pch = 19, 
     xlab = "Days", ylab = "log10(MDS1+10)",
     main = "AMR MDS1")

plot(metaAMR[rownames(MDS), 2], MDS[ ,2], 
     col = "coral3", pch = 19, 
     xlab = "Days", ylab = "log10(MDS2+10)",
     main = "AMR MDS2")

plot(metaAMR[rownames(MDS), 2], MDS[ ,3], 
     col = "coral3", pch = 19, 
     xlab = "Days", ylab = "log10(MDS3+10)",
     main = "AMR MDS3")

plot(metaAMR[rownames(MDS), 2], MDS[ ,4], 
     col = "coral3", pch = 19, 
     xlab = "Days", ylab = "log10(MDS4+10)",
     main = "AMR MDS4")
dev.off()

#AMR
pdf("Plots/TimeVsMDS1234(RGI).pdf", width=12, height=12)
par(mfrow=c(2,2))
par(mar=c(5,6,4,1)+.1)
MDS<-summary(capscale(t(rgiT)~1,distance = "bray"))[["sites"]][,c(1,2,3,4)]
MDS<-log10(MDS+10)
plot(metaRGI[rownames(MDS), 2], MDS[ ,1], 
     col = "cornflowerblue", pch = 19, 
     xlab = "Days", ylab = "log10(MDS1+10)",
     main = "RGI MDS1")

plot(metaRGI[rownames(MDS), 2], MDS[ ,2], 
     col = "cornflowerblue", pch = 19, 
     xlab = "Days", ylab = "log10(MDS2+10)",
     main = "RGI MDS2")

plot(metaRGI[rownames(MDS), 2], MDS[ ,3], 
     col = "cornflowerblue", pch = 19, 
     xlab = "Days", ylab = "log10(MDS3+10)",
     main = "RGI MDS3")

plot(metaRGI[rownames(MDS), 2], MDS[ ,4], 
     col = "cornflowerblue", pch = 19, 
     xlab = "Days", ylab = "log10(MDS4+10)",
     main = "RGI MDS4")
dev.off()

#vsearch
pdf("Plots/TimeVsMDS1234(vsearch).pdf", width=12, height=12)
par(mfrow=c(2,2))
par(mar=c(5,6,4,1)+.1)
MDS<-summary(capscale(t(vsearchT)~1,distance = "bray"))[["sites"]][,c(1,2,3,4)]
MDS<-log10(MDS+10)
plot(metaVSEARCH[rownames(MDS), 2], MDS[ ,1], 
     col = "olivedrab4", pch = 19, 
     xlab = "Days", ylab = "log10(MDS1+10)",
     main = "vsearch MDS1")

plot(metaVSEARCH[rownames(MDS), 2], MDS[ ,2], 
     col = "olivedrab4", pch = 19, 
     xlab = "Days", ylab = "log10(MDS2+10)",
     main = "vsearch MDS2")

plot(metaVSEARCH[rownames(MDS), 2], MDS[ ,3], 
     col = "olivedrab4", pch = 19, 
     xlab = "Days", ylab = "log10(MDS3+10)",
     main = "vsearch MDS3")

plot(metaVSEARCH[rownames(MDS), 2], MDS[ ,4], 
     col = "olivedrab4", pch = 19, 
     xlab = "Days", ylab = "log10(MDS4+10)",
     main = "vsearch MDS4")
dev.off()