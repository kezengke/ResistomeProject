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

#mixed lm 2nd order
pdf("Plots/SortedMixedLmForInOutPatientsTwoColors(Bracken).pdf", width=12, height=18)
par(mfrow=c(3, 2))
par(mar=c(5, 6, 4, 1)+.1)
pvals<-c()
for (i in 1:nrow(brackenT)) {
  myM<-data.frame(unlist(brackenT[i, ]), metaBRACKEN$Timepoint, metaBRACKEN$ptInOut, metaBRACKEN$ID)
  colnames(myM)<-c("counts", "timePoint", "InOut", "ID")
  Model<-lme(counts ~ timePoint * InOut, random = ~1 | ID, data = myM)
  pvals<-rbind(pvals, summary(Model)$tTable[c(2:4), 5])
}
pvals<-data.frame(pvals)
pvals$adjtimePoint<-p.adjust(pvals$timePoint, method = "BH")
pvals$adjInOut<-p.adjust(pvals$InOutOutpatient, method = "BH")
pvals$adjtimePointInOut<-p.adjust(pvals$timePoint.InOutOutpatient, method = "BH")

plot_order<-order(pvals$timePoint.InOutOutpatient)
for (i in 1:nrow(brackenT)) {
  mainText<-paste0(rownames(brackenT)[plot_order[i]], 
                   "\nTimePoint P=", signif(pvals[plot_order[i], 1], 3), ", ", "adj. TimePoint P=", signif(pvals[plot_order[i], 4], 3),
                   "\nInOut P=", signif(pvals[plot_order[i], 2], 3), ", ", "adj. InOut P=", signif(pvals[plot_order[i], 5], 3),
                   "\nTimePoint:InOut P=", signif(pvals[plot_order[i], 3], 3), ", ", "adj. TimePoint:InOut P=", signif(pvals[plot_order[i], 6], 3))
  myM<-data.frame(unlist(brackenT[plot_order[i], ]), metaBRACKEN$Timepoint, metaBRACKEN$ptInOut, metaBRACKEN$ID)
  colnames(myM)<-c("counts", "timePoint", "InOut", "ID")
  Model<-lme(counts ~ timePoint * InOut, random = ~1 | ID, data = myM)
  ModelIn<-lm(counts ~ timePoint, data = myM[myM$InOut == "Inpatient", , drop = F], x=T)
  ModelOut<-lm(counts ~ timePoint, data = myM[myM$InOut == "Outpatient", , drop = F], x=T)
  
  myColor<-vector()
  for (j in 1:length(myM$InOut)) {
    if (myM$InOut[j] == "Inpatient") {
      myColor[j] = "darkorange"
    }
    if (myM$InOut[j] == "Outpatient") {
      myColor[j] = "cornflowerblue"
    }
  }
  plot(metaBRACKEN$Timepoint, unlist(brackenT[plot_order[i],]), 
       main = mainText, cex.main = 1, pch = 19, col = myColor,
       xlab = "TimePoints", ylab = "Counts(Log10)")

  lines(sort(metaBRACKEN$Timepoint), fitted(Model)[order(metaBRACKEN$Timepoint)], 
        col = "grey", type = "l")
  abline(ModelIn, col = "darkorange")
  abline(ModelOut, col = "cornflowerblue")
  
  legend("bottomright",                                     
         legend = c("Mixed Model", "Inpatient", "Outpatient"), 
         lty = 1,                                         
         col = c("grey", "darkorange", "cornflowerblue"), 
         cex = 0.8)                                       

}
dev.off()

pdf("Plots/SortedMixedLmForInOutPatientsTwoColors(AMR).pdf", width=12, height=18)
par(mfrow=c(3, 2))
par(mar=c(5, 6, 4, 1)+.1)
pvals<-c()
for (i in 1:nrow(amrT)) {
  myM<-data.frame(unlist(amrT[i, ]), metaAMR$Timepoint, metaAMR$ptInOut, metaAMR$ID)
  colnames(myM)<-c("counts", "timePoint", "InOut", "ID")
  Model<-lme(counts ~ timePoint * InOut, random = ~1 | ID, data = myM)
  pvals<-rbind(pvals, summary(Model)$tTable[c(2:4), 5])
}
pvals<-data.frame(pvals)
pvals$adjtimePoint<-p.adjust(pvals$timePoint, method = "BH")
pvals$adjInOut<-p.adjust(pvals$InOutOutpatient, method = "BH")
pvals$adjtimePointInOut<-p.adjust(pvals$timePoint.InOutOutpatient, method = "BH")

plot_order<-order(pvals$timePoint.InOutOutpatient)
for (i in 1:nrow(amrT)) {
  mainText<-paste0(rownames(amrT)[plot_order[i]], 
                   "\nTimePoint P=", signif(pvals[plot_order[i], 1], 3), ", ", "adj. TimePoint P=", signif(pvals[plot_order[i], 4], 3),
                   "\nInOut P=", signif(pvals[plot_order[i], 2], 3), ", ", "adj. InOut P=", signif(pvals[plot_order[i], 5], 3),
                   "\nTimePoint:InOut P=", signif(pvals[plot_order[i], 3], 3), ", ", "adj. TimePoint:InOut P=", signif(pvals[plot_order[i], 6], 3))
  myM<-data.frame(unlist(amrT[plot_order[i], ]), metaAMR$Timepoint, metaAMR$ptInOut, metaAMR$ID)
  colnames(myM)<-c("counts", "timePoint", "InOut", "ID")
  Model<-lme(counts ~ timePoint * InOut, random = ~1 | ID, data = myM)
  ModelIn<-lm(counts ~ timePoint, data = myM[myM$InOut == "Inpatient", , drop = F], x=T)
  ModelOut<-lm(counts ~ timePoint, data = myM[myM$InOut == "Outpatient", , drop = F], x=T)
  
  myColor<-vector()
  for (j in 1:length(myM$InOut)) {
    if (myM$InOut[j] == "Inpatient") {
      myColor[j] = "darkorange"
    }
    if (myM$InOut[j] == "Outpatient") {
      myColor[j] = "cornflowerblue"
    }
  }
  plot(metaAMR$Timepoint, unlist(amrT[plot_order[i],]), 
       main = mainText, cex.main = 1, pch = 19, col = myColor,
       xlab = "TimePoints", ylab = "Counts(Log10)")
  
  lines(sort(metaAMR$Timepoint), fitted(Model)[order(metaAMR$Timepoint)], 
        col = "grey", type = "l")
  abline(ModelIn, col = "darkorange")
  abline(ModelOut, col = "cornflowerblue")
  
  legend("bottomright",                                     
         legend = c("Mixed Model", "Inpatient", "Outpatient"), 
         lty = 1,                                         
         col = c("grey", "darkorange", "cornflowerblue"), 
         cex = 0.8)  
}
dev.off()

pdf("Plots/SortedMixedLmForInOutPatientsTwoColors(RGI).pdf", width=12, height=18)
par(mfrow=c(3, 2))
par(mar=c(5, 6, 4, 1)+.1)
pvals<-c()
for (i in 1:nrow(rgiT)) {
  myM<-data.frame(unlist(rgiT[i, ]), metaRGI$Timepoint, metaRGI$ptInOut, metaRGI$ID)
  colnames(myM)<-c("counts", "timePoint", "InOut", "ID")
  Model<-lme(counts ~ timePoint * InOut, random = ~1 | ID, data = myM)
  pvals<-rbind(pvals, summary(Model)$tTable[c(2:4), 5])
}
pvals<-data.frame(pvals)
pvals$adjtimePoint<-p.adjust(pvals$timePoint, method = "BH")
pvals$adjInOut<-p.adjust(pvals$InOutOutpatient, method = "BH")
pvals$adjtimePointInOut<-p.adjust(pvals$timePoint.InOutOutpatient, method = "BH")

plot_order<-order(pvals$timePoint.InOutOutpatient)
for (i in 1:nrow(rgiT)) {
  mainText<-paste0(rownames(rgiT)[plot_order[i]], 
                   "\nTimePoint P=", signif(pvals[plot_order[i], 1], 3), ", ", "adj. TimePoint P=", signif(pvals[plot_order[i], 4], 3),
                   "\nInOut P=", signif(pvals[plot_order[i], 2], 3), ", ", "adj. InOut P=", signif(pvals[plot_order[i], 5], 3),
                   "\nTimePoint:InOut P=", signif(pvals[plot_order[i], 3], 3), ", ", "adj. TimePoint:InOut P=", signif(pvals[plot_order[i], 6], 3))
  myM<-data.frame(unlist(rgiT[plot_order[i], ]), metaRGI$Timepoint, metaRGI$ptInOut, metaRGI$ID)
  colnames(myM)<-c("counts", "timePoint", "InOut", "ID")
  Model<-lme(counts ~ timePoint * InOut, random = ~1 | ID, data = myM)
  ModelIn<-lm(counts ~ timePoint, data = myM[myM$InOut == "Inpatient", , drop = F], x=T)
  ModelOut<-lm(counts ~ timePoint, data = myM[myM$InOut == "Outpatient", , drop = F], x=T)
  
  myColor<-vector()
  for (j in 1:length(myM$InOut)) {
    if (myM$InOut[j] == "Inpatient") {
      myColor[j] = "darkorange"
    }
    if (myM$InOut[j] == "Outpatient") {
      myColor[j] = "cornflowerblue"
    }
  }
  plot(metaRGI$Timepoint, unlist(rgiT[plot_order[i],]), 
       main = mainText, cex.main = 1, pch = 19, col = myColor,
       xlab = "TimePoints", ylab = "Counts(Log10)")
  
  lines(sort(metaRGI$Timepoint), fitted(Model)[order(metaRGI$Timepoint)], 
        col = "grey", type = "l")
  abline(ModelIn, col = "darkorange")
  abline(ModelOut, col = "cornflowerblue")
  
  legend("bottomright",                                     
         legend = c("Mixed Model", "Inpatient", "Outpatient"), 
         lty = 1,                                         
         col = c("grey", "darkorange", "cornflowerblue"), 
         cex = 0.8)  
}
dev.off()

pdf("Plots/SortedMixedLmForInOutPatientsTwoColors(vsearch).pdf", width=12, height=18)
par(mfrow=c(3, 2))
par(mar=c(5, 6, 4, 1)+.1)
pvals<-c()
for (i in 1:nrow(vsearchT)) {
  myM<-data.frame(unlist(vsearchT[i, ]), metaVSEARCH$Timepoint, metaVSEARCH$ptInOut, metaVSEARCH$ID)
  colnames(myM)<-c("counts", "timePoint", "InOut", "ID")
  Model<-lme(counts ~ timePoint * InOut, random = ~1 | ID, data = myM)
  pvals<-rbind(pvals, summary(Model)$tTable[c(2:4), 5])
}
pvals<-data.frame(pvals)
pvals$adjtimePoint<-p.adjust(pvals$timePoint, method = "BH")
pvals$adjInOut<-p.adjust(pvals$InOutOutpatient, method = "BH")
pvals$adjtimePointInOut<-p.adjust(pvals$timePoint.InOutOutpatient, method = "BH")

plot_order<-order(pvals$timePoint.InOutOutpatient)
for (i in 1:nrow(vsearchT)) {
  mainText<-paste0(rownames(vsearchT)[plot_order[i]], 
                   "\nTimePoint P=", signif(pvals[plot_order[i], 1], 3), ", ", "adj. TimePoint P=", signif(pvals[plot_order[i], 4], 3),
                   "\nInOut P=", signif(pvals[plot_order[i], 2], 3), ", ", "adj. InOut P=", signif(pvals[plot_order[i], 5], 3),
                   "\nTimePoint:InOut P=", signif(pvals[plot_order[i], 3], 3), ", ", "adj. TimePoint:InOut P=", signif(pvals[plot_order[i], 6], 3))
  myM<-data.frame(unlist(vsearchT[plot_order[i], ]), metaVSEARCH$Timepoint, metaVSEARCH$ptInOut, metaVSEARCH$ID)
  colnames(myM)<-c("counts", "timePoint", "InOut", "ID")
  Model<-lme(counts ~ timePoint * InOut, random = ~1 | ID, data = myM)
  ModelIn<-lm(counts ~ timePoint, data = myM[myM$InOut == "Inpatient", , drop = F], x=T)
  ModelOut<-lm(counts ~ timePoint, data = myM[myM$InOut == "Outpatient", , drop = F], x=T)
  
  myColor<-vector()
  for (j in 1:length(myM$InOut)) {
    if (myM$InOut[j] == "Inpatient") {
      myColor[j] = "darkorange"
    }
    if (myM$InOut[j] == "Outpatient") {
      myColor[j] = "cornflowerblue"
    }
  }
  plot(metaVSEARCH$Timepoint, unlist(vsearchT[plot_order[i],]), 
       main = mainText, cex.main = 1, pch = 19, col = myColor,
       xlab = "TimePoints", ylab = "Counts(Log10)")
  
  lines(sort(metaVSEARCH$Timepoint), fitted(Model)[order(metaVSEARCH$Timepoint)], 
        col = "grey", type = "l")
  abline(ModelIn, col = "darkorange")
  abline(ModelOut, col = "cornflowerblue")
  
  legend("bottomright",                                     
         legend = c("Mixed Model", "Inpatient", "Outpatient"), 
         lty = 1,                                         
         col = c("grey", "darkorange", "cornflowerblue"), 
         cex = 0.8)  
}
dev.off()
