rm(list = ls())
library("nlme")

metaData<-read.csv("metaWithBins.csv", header = T, row.names = 1)

#counts tables
brackenT<-read.csv("CountsTables/brackenNormalized.csv", header = T, row.names = 1, check.names = F)
amrT<-read.csv("CountsTables/amrNormalized.csv", header = T, row.names = 1, check.names = F)
rgiT<-read.csv("CountsTables/rgiNormalized.csv", header = T, row.names = 1, check.names = F)
vsearchT<-read.csv("CountsTables/vsearchNormalized.csv", header = T, row.names = 1, check.names = F)

#filtering
brackenT<-brackenT[apply(brackenT == 0, 1, sum) <= (ncol(brackenT)*0.8), ]
amrT<-amrT[apply(amrT == 0, 1, sum) <= (ncol(amrT)*0.8), ]
rgiT<-rgiT[apply(rgiT == 0, 1, sum) <= (ncol(rgiT)*0.8), ]
vsearchT<-vsearchT[apply(vsearchT == 0, 1, sum) <= (ncol(vsearchT)*0.8), ]

metaBRACKEN<-metaData[colnames(brackenT), , drop = F]
metaAMR<-metaData[colnames(amrT), , drop = F]
metaRGI<-metaData[colnames(rgiT), , drop = F]
metaVSEARCH<-metaData[colnames(vsearchT), , drop = F]

#mixed lm 2nd order
pdf("Plots/SortedSecondOrderMixedLm(Bracken).pdf", width=12, height=18)
par(mfrow=c(3, 2))
par(mar=c(5, 6, 4, 1)+.1)
pvals<-vector()
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
  pvals[i]<-anova(polyM)[2,4]
}
adjPvals<-p.adjust(pvals, method = "BH")
adjFpvals<-p.adjust(Fpval, method = "BH")

plot_order<-order(pvals)
for (i in 1:nrow(brackenT)) {
  mainText<-paste0(rownames(brackenT)[plot_order[i]], 
                   "\nANOVA P=", signif(pvals[plot_order[i]], 3), ", ", "adj. ANOVA-P=", signif(adjPvals[plot_order[i]], 3), 
                   "\nLine vs. Poly F-P=", signif(Fpval[plot_order[i]], 3), ", ", "adj. F-P=", signif(adjFpvals[plot_order[i]], 3))
  myM<-data.frame(unlist(brackenT[plot_order[i], ]), metaBRACKEN$Timepoint, metaBRACKEN$ID)
  colnames(myM)<-c("counts", "timePoint", "ID")
  Model<-lme(counts ~ poly(timePoint, 2), random = ~1 | ID, data = myM)
  plot(metaBRACKEN$Timepoint, unlist(brackenT[plot_order[i],]), 
       pch = 19, col = "tan2", 
       main = mainText,
       xlab = "TimePoints", ylab = "Counts(Log10)")
  lines(sort(metaBRACKEN$Timepoint), fitted(Model)[order(metaBRACKEN$Timepoint)], 
        col = "grey", type = "l")
}
dev.off()

pdf("Plots/SortedSecondOrderMixedLm(AMR).pdf", width=12, height=18)
par(mfrow=c(3, 2))
par(mar=c(5, 6, 4, 1)+.1)
pvals<-vector()
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
  pvals[i]<-anova(polyM)[2,4]
}
adjPvals<-p.adjust(pvals, method = "BH")
adjFpvals<-p.adjust(Fpval, method = "BH")

plot_order<-order(pvals)
for (i in 1:nrow(amrT)) {
  mainText<-paste0(rownames(amrT)[plot_order[i]], "\nANOVA P=", signif(pvals[plot_order[i]], 3), ", ", "adj. ANOVA-P=", signif(adjPvals[plot_order[i]], 3), 
                   "\nLine vs. Poly F-P=", signif(Fpval[plot_order[i]], 3), ", ", "adj. F-P=", signif(adjFpvals[plot_order[i]], 3))
  myM<-data.frame(unlist(amrT[plot_order[i], ]), metaAMR$Timepoint, metaAMR$ID)
  colnames(myM)<-c("counts", "timePoint", "ID")
  Model<-lme(counts ~ poly(timePoint, 2), random = ~1 | ID, data = myM)
  plot(metaAMR$Timepoint, unlist(amrT[plot_order[i],]), 
       pch = 19, col = "coral3", 
       main = mainText,
       xlab = "TimePoints", ylab = "Counts(Log10)")
  lines(sort(metaAMR$Timepoint), fitted(Model)[order(metaAMR$Timepoint)], 
        col = "grey", type = "l")
}
dev.off()

pdf("Plots/SortedSecondOrderMixedLm(RGI).pdf", width=12, height=18)
par(mfrow=c(3, 2))
par(mar=c(5, 6, 4, 1)+.1)
rgiT<-rgiT[apply(rgiT == 0, 1, sum) <= (ncol(rgiT)*0.8), ]
pvals<-vector()
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
  pvals[i]<-anova(polyM)[2,4]
}
adjPvals<-p.adjust(pvals, method = "BH")
adjFpvals<-p.adjust(Fpval, method = "BH")

plot_order<-order(pvals)
for (i in 1:nrow(rgiT)) {
  mainText<-paste0(rownames(rgiT)[plot_order[i]], "\nANOVA P=", signif(pvals[plot_order[i]], 3), ", ", "adj. ANOVA-P=", signif(adjPvals[plot_order[i]], 3), 
                   "\nLine vs. Poly F-P=", signif(Fpval[plot_order[i]], 3), ", ", "adj. F-P=", signif(adjFpvals[plot_order[i]], 3))
  myM<-data.frame(unlist(rgiT[plot_order[i], ]), metaRGI$Timepoint, metaRGI$ID)
  colnames(myM)<-c("counts", "timePoint", "ID")
  Model<-lme(counts ~ poly(timePoint, 2), random = ~1 | ID, data = myM)
  plot(metaRGI$Timepoint, unlist(rgiT[plot_order[i],]), 
       pch = 19, col = "cornflowerblue", 
       main = mainText,
       xlab = "TimePoints", ylab = "Counts(Log10)")
  lines(sort(metaRGI$Timepoint), fitted(Model)[order(metaRGI$Timepoint)], 
        col = "grey", type = "l")
}
dev.off()

pdf("Plots/SortedSecondOrderMixedLm(vsearch).pdf", width=12, height=18)
par(mfrow=c(3, 2))
par(mar=c(5, 6, 4, 1)+.1)
pvals<-vector()
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
  pvals[i]<-anova(polyM)[2,4]
}
adjPvals<-p.adjust(pvals, method = "BH")
adjFpvals<-p.adjust(Fpval, method = "BH")

plot_order<-order(pvals)
for (i in 1:nrow(vsearchT)) {
  mainText<-paste0(strsplit(rownames(vsearchT)[plot_order[i]], "\\|")[[1]][6],  "\nANOVA P=", signif(pvals[plot_order[i]], 3), ", ", "adj. ANOVA-P=", signif(adjPvals[plot_order[i]], 3), 
                   "\nLine vs. Poly F-P=", signif(Fpval[plot_order[i]], 3), ", ", "adj. F-P=", signif(adjFpvals[plot_order[i]], 3))
  myM<-data.frame(unlist(vsearchT[plot_order[i], ]), metaVSEARCH$Timepoint, metaVSEARCH$ID)
  colnames(myM)<-c("counts", "timePoint", "ID")
  Model<-lme(counts ~ poly(timePoint, 2), random = ~1 | ID, data = myM)
  plot(metaVSEARCH$Timepoint, unlist(vsearchT[plot_order[i],]), 
       pch = 19, col = "olivedrab4", 
       main = mainText,
       xlab = "TimePoints", ylab = "Counts(Log10)")
  lines(sort(metaVSEARCH$Timepoint), fitted(Model)[order(metaVSEARCH$Timepoint)], 
        col = "grey", type = "l")
}
dev.off()
