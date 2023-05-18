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

pdf("Plots/MixedLmForPatientIDPvalHistograms.pdf", width=12, height=12)
par(mfrow=c(2, 2))
par(mar=c(5, 6, 4, 1)+.1)

modelPvals<-vector()
for (i in 1:nrow(brackenT)) {
  myM<-data.frame(unlist(brackenT[i, ]), metaBRACKEN$Timepoint, metaBRACKEN$ID)
  colnames(myM)<-c("counts", "timePoint", "ID")
  Model<-lme(counts ~ timePoint, random = ~1 | ID, data = myM)
  modelPvals[i]<-anova(Model)[2,4]
  
}
hist(modelPvals, breaks=seq(0, 1, 0.05), xlab = "p-value", col = "tan2",
     main = "Bracken(Species) Mixed LM ANOVA P-vals", 
     xlim=c(0,1), cex.lab = 1.5, cex.main = 1.7, cex.axis = 1.4)
abline(v=0.05, col=gray(.5), lty=2)
mtext(at=0.05, side=3, text=0.05, col=gray(.5))

modelPvals<-vector()
for (i in 1:nrow(amrT)) {
  myM<-data.frame(unlist(amrT[i, ]), metaAMR$Timepoint, metaAMR$ID)
  colnames(myM)<-c("counts", "timePoint", "ID")
  Model<-lme(counts ~ timePoint, random = ~1 | ID, data = myM)
  modelPvals[i]<-anova(Model)[2,4]
  
}
hist(modelPvals, breaks=seq(0, 1, 0.05), xlab = "p-value", col = "coral3",
     main = "AMR Mixed LM ANOVA P-vals", 
     xlim=c(0,1), cex.lab = 1.5, cex.main = 1.7, cex.axis = 1.4)
abline(v=0.05, col=gray(.5), lty=2)
mtext(at=0.05, side=3, text=0.05, col=gray(.5))

modelPvals<-vector()
for (i in 1:nrow(rgiT)) {
  myM<-data.frame(unlist(rgiT[i, ]), metaRGI$Timepoint, metaRGI$ID)
  colnames(myM)<-c("counts", "timePoint", "ID")
  Model<-lme(counts ~ timePoint, random = ~1 | ID, data = myM)
  modelPvals[i]<-anova(Model)[2,4]
  
}
hist(modelPvals, breaks=seq(0, 1, 0.05), xlab = "p-value", col = "cornflowerblue",
     main = "RGI Mixed LM ANOVA P-vals", 
     xlim=c(0,1), cex.lab = 1.5, cex.main = 1.7, cex.axis = 1.4)
abline(v=0.05, col=gray(.5), lty=2)
mtext(at=0.05, side=3, text=0.05, col=gray(.5))

modelPvals<-vector()
for (i in 1:nrow(vsearchT)) {
  myM<-data.frame(unlist(vsearchT[i, ]), metaVSEARCH$Timepoint, metaVSEARCH$ID)
  colnames(myM)<-c("counts", "timePoint", "ID")
  Model<-lme(counts ~ timePoint, random = ~1 | ID, data = myM)
  modelPvals[i]<-anova(Model)[2,4]
  
}
hist(modelPvals, breaks=seq(0, 1, 0.05), xlab = "p-value", col = "olivedrab4",
     main = "vsearch Mixed LM ANOVA P-vals", 
     xlim=c(0,1), cex.lab = 1.5, cex.main = 1.7, cex.axis = 1.4)
abline(v=0.05, col=gray(.5), lty=2)
mtext(at=0.05, side=3, text=0.05, col=gray(.5))

dev.off()
