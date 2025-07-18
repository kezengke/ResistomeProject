rm(list = ls())
library("nlme")

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

pdf("Plots/LineVsPolyFstatsPvalHistograms.pdf", width=12, height=12)
par(mfrow=c(2, 2))
par(mar=c(5, 6, 4, 1)+.1)

Fpval<-vector()
for (i in 1:nrow(brackenT)) {
  myM<-data.frame(unlist(brackenT[i, ]), metaBRACKEN$Timepoint, metaBRACKEN$ID)
  colnames(myM)<-c("counts", "timePoint", "ID")
  
  lineM<-lme(counts ~ timePoint, random = ~1 | ID, data = myM)
  polyM<-lme(counts ~ poly(timePoint, 2), random = ~1 | ID, data = myM)
  
  reducedError<-sum(resid(lineM)^2)
  fullError<-sum(resid(polyM)^2)
  
  reducedDF<-nrow(myM)-2
  fullDF<-nrow(myM)-3
  
  myF <- ((reducedError - fullError)/(reducedDF - fullDF))/(fullError/fullDF)
  
  Fpval[i]<-pf(myF, 1, fullDF, lower.tail = F)
}

hist(Fpval, breaks=seq(0, 1, 0.05), xlab = "p-value", col = "tan2",
     main = paste("Bracken(Species)", "\nline vs. poly F-test P-vals"), 
     xlim=c(0,1), cex.lab = 1.5, cex.main = 1.7, cex.axis = 1.4)
abline(v=0.05, col=gray(.5), lty=2)
mtext(at=0.05, side=3, text=0.05, col=gray(.5))

Fpval<-vector()
for (i in 1:nrow(amrT)) {
  myM<-data.frame(unlist(amrT[i, ]), metaAMR$Timepoint, metaAMR$ID)
  colnames(myM)<-c("counts", "timePoint", "ID")
  
  lineM<-lme(counts ~ timePoint, random = ~1 | ID, data = myM)
  polyM<-lme(counts ~ poly(timePoint, 2), random = ~1 | ID, data = myM)
  
  reducedError<-sum(resid(lineM)^2)
  fullError<-sum(resid(polyM)^2)
  
  reducedDF<-nrow(myM)-2
  fullDF<-nrow(myM)-3
  
  myF <- ((reducedError - fullError)/(reducedDF - fullDF))/(fullError/fullDF)
  
  Fpval[i]<-pf(myF, 1, fullDF, lower.tail = F)
}

hist(Fpval, breaks=seq(0, 1, 0.05), xlab = "p-value", col = "coral3",
     main = paste("AMR", "\nline vs. poly F-test P-vals"), 
     xlim=c(0,1), cex.lab = 1.5, cex.main = 1.7, cex.axis = 1.4)
abline(v=0.05, col=gray(.5), lty=2)
mtext(at=0.05, side=3, text=0.05, col=gray(.5))

Fpval<-vector()
for (i in 1:nrow(rgiT)) {
  myM<-data.frame(unlist(rgiT[i, ]), metaRGI$Timepoint, metaRGI$ID)
  colnames(myM)<-c("counts", "timePoint", "ID")
  
  lineM<-lme(counts ~ timePoint, random = ~1 | ID, data = myM)
  polyM<-lme(counts ~ poly(timePoint, 2), random = ~1 | ID, data = myM)
  
  reducedError<-sum(resid(lineM)^2)
  fullError<-sum(resid(polyM)^2)
  
  reducedDF<-nrow(myM)-2
  fullDF<-nrow(myM)-3
  
  myF <- ((reducedError - fullError)/(reducedDF - fullDF))/(fullError/fullDF)
  
  Fpval[i]<-pf(myF, 1, fullDF, lower.tail = F)
}

hist(Fpval, breaks=seq(0, 1, 0.05), xlab = "p-value", col = "cornflowerblue",
     main = paste("RGI", "\nline vs. poly F-test P-vals"), 
     xlim=c(0,1), cex.lab = 1.5, cex.main = 1.7, cex.axis = 1.4)
abline(v=0.05, col=gray(.5), lty=2)
mtext(at=0.05, side=3, text=0.05, col=gray(.5))

Fpval<-vector()
for (i in 1:nrow(vsearchT)) {
  myM<-data.frame(unlist(vsearchT[i, ]), metaVSEARCH$Timepoint, metaVSEARCH$ID)
  colnames(myM)<-c("counts", "timePoint", "ID")
  
  lineM<-lme(counts ~ timePoint, random = ~1 | ID, data = myM)
  polyM<-lme(counts ~ poly(timePoint, 2), random = ~1 | ID, data = myM)
  
  reducedError<-sum(resid(lineM)^2)
  fullError<-sum(resid(polyM)^2)
  
  reducedDF<-nrow(myM)-2
  fullDF<-nrow(myM)-3
  
  myF <- ((reducedError - fullError)/(reducedDF - fullDF))/(fullError/fullDF)
  
  Fpval[i]<-pf(myF, 1, fullDF, lower.tail = F)
}

hist(Fpval, breaks=seq(0, 1, 0.05), xlab = "p-value", col = "olivedrab4",
     main = paste("vsearch", "\nline vs. poly F-test P-vals"), 
     xlim=c(0,1), cex.lab = 1.5, cex.main = 1.7, cex.axis = 1.4)
abline(v=0.05, col=gray(.5), lty=2)
mtext(at=0.05, side=3, text=0.05, col=gray(.5))

dev.off()
