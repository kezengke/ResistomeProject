# pval histogram of mixed linear models
rm(list = ls())
library(stringr)
library("dplyr")
library("nlme")

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
amrMatchT<-amrMatchT[, c(18,11)] #keep only the gene name and Corrected AMR classification column

amrTypes<-unique(amrMatchT$AMR.Gene.Family)
newamrT<-c()
for (i in 1:length(amrTypes)) {
  genes<-amrMatchT %>% filter(amrMatchT$AMR.Gene.Family == AMR.Gene.Family[i])
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
rgiMatchT<-rgiMatchT[, c(15,11)] #keep only the gene name and Corrected AMR classification column

amrTypes<-unique(rgiMatchT$AMR.Gene.Family)
newrgiT<-c()
for (i in 1:length(amrTypes)) {
  genes<-rgiMatchT %>% filter(rgiMatchT$AMR.Gene.Family == amrTypes[i])
  currentFrame<-rgiT[genes$RGI.CARD.Short.Name, ]
  sumRow<-colSums(currentFrame)
  newrgiT<-rbind(newrgiT, sumRow)
}
rownames(newrgiT)<-amrTypes

#VSEARCH
#keep only the rows present in vsearch counts table in master list
vsearchMatchT<-geneMatchT %>% filter(row_number() %in% match(rownames(vsearchT), geneMatchT$Vsearch.ARO.Name))
vsearchMatchT<-vsearchMatchT[, c(7,11)] #keep only the gene name and Corrected AMR classification column

amrTypes<-unique(vsearchMatchT$AMR.Gene.Family)
newvsearchT<-c()
for (i in 1:length(amrTypes)) {
  genes<-vsearchMatchT %>% filter(vsearchMatchT$AMR.Gene.Family == amrTypes[i])
  currentFrame<-vsearchT[genes$Vsearch.ARO.Name, ]
  sumRow<-colSums(currentFrame)
  newvsearchT<-rbind(newvsearchT, sumRow)
}
rownames(newvsearchT)<-amrTypes

#filtering
amrT<-amrT[apply(amrT == 0, 1, sum) <= (ncol(amrT)*0.8), ]
rgiT<-rgiT[apply(rgiT == 0, 1, sum) <= (ncol(rgiT)*0.8), ]
vsearchT<-vsearchT[apply(vsearchT == 0, 1, sum) <= (ncol(vsearchT)*0.8), ]

metaAMR<-metaData[colnames(newamrT), , drop = F]
metaRGI<-metaData[colnames(newrgiT), , drop = F]
metaVSEARCH<-metaData[colnames(newvsearchT), , drop = F]

pdf("Plots/MixedSecondOrderLmAMRGeneFamilyPvalHistograms.pdf", width=18, height=6)
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
     main = "AMR(AMRGeneFamily) Mixed LM 2nd Order ANOVA P-vals", 
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
     main = "RGI(AMRGeneFamily) Mixed LM 2nd Order ANOVA P-vals", 
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
     main = "vsearch(AMRGeneFamily) Mixed LM 2nd Order ANOVA P-vals", 
     xlim=c(0,1), cex.lab = 1.5, cex.main = 1.7, cex.axis = 1.4)
abline(v=0.05, col=gray(.5), lty=2)
mtext(at=0.05, side=3, text=0.05, col=gray(.5))

dev.off()
