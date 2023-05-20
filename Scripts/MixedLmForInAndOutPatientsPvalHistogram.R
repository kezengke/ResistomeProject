rm(list = ls())
library("nlme")

metaData<-read.csv("metaWithBins.csv", header = T, row.names = 1)

#counts tables
brackenT<-read.csv("CountsTables/brackenFiltered.csv", header = T, row.names = 1, check.names = F)
amrT<-read.csv("CountsTables/amrFiltered.csv", header = T, row.names = 1, check.names = F)
rgiT<-read.csv("CountsTables/rgiFiltered.csv", header = T, row.names = 1, check.names = F)
vsearchT<-read.csv("CountsTables/vsearchFiltered.csv", header = T, row.names = 1, check.names = F)

#meta for each table
metaBRACKEN<-metaData[colnames(brackenT), , drop = F]
metaAMR<-metaData[colnames(amrT), , drop = F]
metaRGI<-metaData[colnames(rgiT), , drop = F]
metaVSEARCH<-metaData[colnames(vsearchT), , drop = F]

pdf("Plots/MixedLmForInAndOutPatientsPvalHistograms.pdf", width=18, height=24)
par(mfrow=c(4, 3))
par(mar=c(5, 6, 4, 1)+.1)

modelPvals<-c()
for (i in 1:nrow(brackenT)) {
  myM<-data.frame(unlist(brackenT[i, ]), metaBRACKEN$Timepoint, metaBRACKEN$ptInOut, metaBRACKEN$ID)
  colnames(myM)<-c("counts", "timePoint", "InOut", "ID")
  Model<-lme(counts ~ timePoint * InOut, random = ~1 | ID, data = myM)
  modelPvals<-rbind(modelPvals, summary(Model)$tTable[c(2:4), 5])
  
}
hist(modelPvals[,1], breaks=seq(0, 1, 0.05), xlab = "p-value", col = "tan2",
     main = "Bracken(Species) Mixed LM TimePoint P-vals", 
     xlim=c(0,1), cex.lab = 1.5, cex.main = 1.7, cex.axis = 1.4)
abline(v=0.05, col=gray(.5), lty=2)
mtext(at=0.05, side=3, text=0.05, col=gray(.5))

hist(modelPvals[,2], breaks=seq(0, 1, 0.05), xlab = "p-value", col = "tan2",
     main = "Bracken(Species) Mixed LM InOutPatient P-vals", 
     xlim=c(0,1), cex.lab = 1.5, cex.main = 1.7, cex.axis = 1.4)
abline(v=0.05, col=gray(.5), lty=2)
mtext(at=0.05, side=3, text=0.05, col=gray(.5))

hist(modelPvals[,3], breaks=seq(0, 1, 0.05), xlab = "p-value", col = "tan2",
     main = "Bracken(Species) Mixed LM TimePoint:InOutPatient P-vals", 
     xlim=c(0,1), cex.lab = 1.5, cex.main = 1.7, cex.axis = 1.4)
abline(v=0.05, col=gray(.5), lty=2)
mtext(at=0.05, side=3, text=0.05, col=gray(.5))

modelPvals<-c()
for (i in 1:nrow(amrT)) {
  myM<-data.frame(unlist(amrT[i, ]), metaAMR$Timepoint, metaAMR$ptInOut, metaAMR$ID)
  colnames(myM)<-c("counts", "timePoint", "InOut", "ID")
  Model<-lme(counts ~ timePoint * InOut, random = ~1 | ID, data = myM)
  modelPvals<-rbind(modelPvals, summary(Model)$tTable[c(2:4), 5])
  
}
hist(modelPvals[,1], breaks=seq(0, 1, 0.05), xlab = "p-value", col = "coral3",
     main = "AMR Mixed LM TimePoint P-vals", 
     xlim=c(0,1), cex.lab = 1.5, cex.main = 1.7, cex.axis = 1.4)
abline(v=0.05, col=gray(.5), lty=2)
mtext(at=0.05, side=3, text=0.05, col=gray(.5))

hist(modelPvals[,2], breaks=seq(0, 1, 0.05), xlab = "p-value", col = "coral3",
     main = "AMR Mixed LM InOutPatient P-vals", 
     xlim=c(0,1), cex.lab = 1.5, cex.main = 1.7, cex.axis = 1.4)
abline(v=0.05, col=gray(.5), lty=2)
mtext(at=0.05, side=3, text=0.05, col=gray(.5))

hist(modelPvals[,3], breaks=seq(0, 1, 0.05), xlab = "p-value", col = "coral3",
     main = "AMR Mixed LM TimePoint:InOutPatient P-vals", 
     xlim=c(0,1), cex.lab = 1.5, cex.main = 1.7, cex.axis = 1.4)
abline(v=0.05, col=gray(.5), lty=2)
mtext(at=0.05, side=3, text=0.05, col=gray(.5))

modelPvals<-c()
for (i in 1:nrow(rgiT)) {
  myM<-data.frame(unlist(rgiT[i, ]), metaRGI$Timepoint, metaRGI$ptInOut, metaRGI$ID)
  colnames(myM)<-c("counts", "timePoint", "InOut", "ID")
  Model<-lme(counts ~ timePoint * InOut, random = ~1 | ID, data = myM)
  modelPvals<-rbind(modelPvals, summary(Model)$tTable[c(2:4), 5])
  
}
hist(modelPvals[,1], breaks=seq(0, 1, 0.05), xlab = "p-value", col = "cornflowerblue",
     main = "RGI Mixed LM TimePoint P-vals", 
     xlim=c(0,1), cex.lab = 1.5, cex.main = 1.7, cex.axis = 1.4)
abline(v=0.05, col=gray(.5), lty=2)
mtext(at=0.05, side=3, text=0.05, col=gray(.5))

hist(modelPvals[,2], breaks=seq(0, 1, 0.05), xlab = "p-value", col = "cornflowerblue",
     main = "RGI Mixed LM InOutPatient P-vals", 
     xlim=c(0,1), cex.lab = 1.5, cex.main = 1.7, cex.axis = 1.4)
abline(v=0.05, col=gray(.5), lty=2)
mtext(at=0.05, side=3, text=0.05, col=gray(.5))

hist(modelPvals[,3], breaks=seq(0, 1, 0.05), xlab = "p-value", col = "cornflowerblue",
     main = "RGI Mixed LM TimePoint:InOutPatient P-vals", 
     xlim=c(0,1), cex.lab = 1.5, cex.main = 1.7, cex.axis = 1.4)
abline(v=0.05, col=gray(.5), lty=2)
mtext(at=0.05, side=3, text=0.05, col=gray(.5))

modelPvals<-c()
for (i in 1:nrow(vsearchT)) {
  myM<-data.frame(unlist(vsearchT[i, ]), metaVSEARCH$Timepoint, metaVSEARCH$ptInOut, metaVSEARCH$ID)
  colnames(myM)<-c("counts", "timePoint", "InOut", "ID")
  Model<-lme(counts ~ timePoint * InOut, random = ~1 | ID, data = myM)
  modelPvals<-rbind(modelPvals, summary(Model)$tTable[c(2:4), 5])
  
}
hist(modelPvals[,1], breaks=seq(0, 1, 0.05), xlab = "p-value", col = "olivedrab4",
     main = "vsearch Mixed LM TimePoint P-vals", 
     xlim=c(0,1), cex.lab = 1.5, cex.main = 1.7, cex.axis = 1.4)
abline(v=0.05, col=gray(.5), lty=2)
mtext(at=0.05, side=3, text=0.05, col=gray(.5))

hist(modelPvals[,2], breaks=seq(0, 1, 0.05), xlab = "p-value", col = "olivedrab4",
     main = "vsearch Mixed LM InOutPatient P-vals", 
     xlim=c(0,1), cex.lab = 1.5, cex.main = 1.7, cex.axis = 1.4)
abline(v=0.05, col=gray(.5), lty=2)
mtext(at=0.05, side=3, text=0.05, col=gray(.5))

hist(modelPvals[,3], breaks=seq(0, 1, 0.05), xlab = "p-value", col = "olivedrab4",
     main = "vsearch Mixed LM TimePoint:InOutPatient P-vals", 
     xlim=c(0,1), cex.lab = 1.5, cex.main = 1.7, cex.axis = 1.4)
abline(v=0.05, col=gray(.5), lty=2)
mtext(at=0.05, side=3, text=0.05, col=gray(.5))

dev.off()
