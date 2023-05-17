# pval histogram of mixed linear models
rm(list = ls())
library(stringr)
library("dplyr")
library("nlme")

geneMatchT<-read.csv("MASTER_AMRlist_2023_03.csv", sep = ",", header = T)
metaData<-read.csv("metaWithBins.csv", header = T, row.names = 1)

#gene counts tables
amrT<-read.csv("CountsTables/amrFiltered.csv", header = T, row.names = 1, check.names = F)
rgiT<-read.csv("CountsTables/rgiFiltered.csv", header = T, row.names = 1, check.names = F)
vsearchT<-read.csv("CountsTables/vsearchFiltered.csv", header = T, row.names = 1, check.names = F)

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

amrTypes<-unique(vsearchMatchT$Corrected.AMR.classification)
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

pdf("Plots/MixedSecondOrderLmCorrectedAMRClassificationPvalHistograms.pdf", width=18, height=6)
par(mfrow=c(1, 3))
par(mar=c(5, 6, 4, 1)+.1)

modelPvals<-vector()
for (i in 1:nrow(newamrT)) {
  myM<-data.frame(unlist(newamrT[i, ]), metaAMR$Timepoint, metaAMR$ID)
  colnames(myM)<-c("counts", "timePoint", "ID")
  Model<-lme(counts ~ poly(timePoint, 2), random = ~1 | ID, data = myM)
  modelPvals[i]<-anova(Model)[2,4]
  
}
hist(modelPvals, breaks=seq(0, 1, 0.05), xlab = "p-value", col = "coral3",
     main = "AMR(CorrectedAMRClassification) Mixed LM 2nd Order ANOVA P-vals", 
     xlim=c(0,1), cex.lab = 1.5, cex.main = 1.7, cex.axis = 1.4)
abline(v=0.05, col=gray(.5), lty=2)
mtext(at=0.05, side=3, text=0.05, col=gray(.5))

modelPvals<-vector()
for (i in 1:nrow(newrgiT)) {
  myM<-data.frame(unlist(newrgiT[i, ]), metaRGI$Timepoint, metaRGI$ID)
  colnames(myM)<-c("counts", "timePoint", "ID")
  Model<-lme(counts ~ poly(timePoint, 2), random = ~1 | ID, data = myM)
  modelPvals[i]<-anova(Model)[2,4]
  
}
hist(modelPvals, breaks=seq(0, 1, 0.05), xlab = "p-value", col = "cornflowerblue",
     main = "RGI(CorrectedAMRClassification) Mixed LM 2nd Order ANOVA P-vals", 
     xlim=c(0,1), cex.lab = 1.5, cex.main = 1.7, cex.axis = 1.4)
abline(v=0.05, col=gray(.5), lty=2)
mtext(at=0.05, side=3, text=0.05, col=gray(.5))

modelPvals<-vector()
for (i in 1:nrow(newvsearchT)) {
  myM<-data.frame(unlist(newvsearchT[i, ]), metaVSEARCH$Timepoint, metaVSEARCH$ID)
  colnames(myM)<-c("counts", "timePoint", "ID")
  Model<-lme(counts ~ poly(timePoint, 2), random = ~1 | ID, data = myM)
  modelPvals[i]<-anova(Model)[2,4]
  
}
hist(modelPvals, breaks=seq(0, 1, 0.05), xlab = "p-value", col = "olivedrab4",
     main = "vsearch(CorrectedAMRClassification) Mixed LM 2nd Order ANOVA P-vals", 
     xlim=c(0,1), cex.lab = 1.5, cex.main = 1.7, cex.axis = 1.4)
abline(v=0.05, col=gray(.5), lty=2)
mtext(at=0.05, side=3, text=0.05, col=gray(.5))

dev.off()
