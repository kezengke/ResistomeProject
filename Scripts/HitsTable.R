rm(list = ls())
library(stringr)
library("dplyr")
library("nlme")

#gene category match table
geneMatchT<-read.csv("MASTER_AMRlist_2023_03.csv", sep = ",", header = T)

#metadata
metaData<-read.csv("metaWithBins.csv", header = T, row.names = 1)

#counts tables
amrT<-read.csv("CountsTables/amrNormalized.csv", header = T, row.names = 1, check.names = F)
rgiT<-read.csv("CountsTables/rgiNormalized.csv", header = T, row.names = 1, check.names = F)
vsearchT<-read.csv("CountsTables/vsearchNormalized.csv", header = T, row.names = 1, check.names = F)

#filtering
amrT<-amrT[apply(amrT == 0, 1, sum) <= (ncol(amrT)*0.8), ]
rgiT<-rgiT[apply(rgiT == 0, 1, sum) <= (ncol(rgiT)*0.8), ]
vsearchT<-vsearchT[apply(vsearchT == 0, 1, sum) <= (ncol(vsearchT)*0.8), ]

#meta for each table
metaAMR<-metaData[colnames(amrT), , drop = F]
metaRGI<-metaData[colnames(rgiT), , drop = F]
metaVSEARCH<-metaData[colnames(vsearchT), , drop = F]

bigT<-c()
#mixed lm 2nd order
amrpval<-vector()
for (i in 1:nrow(amrT)) {
  myM<-data.frame(unlist(amrT[i, ]), metaAMR$Timepoint, metaAMR$ID)
  colnames(myM)<-c("counts", "timePoint", "ID")
  Model<-lme(counts ~ poly(timePoint, 2), random = ~1 | ID, data = myM)
  amrpval[i]<-anova(Model)[2,4]
}
adjamrpval<-p.adjust(amrpval, method = "BH")
amr<-c(sum(amrpval<0.05), sum(adjamrpval<0.05))

rgipval<-vector()
for (i in 1:nrow(rgiT)) {
  myM<-data.frame(unlist(rgiT[i, ]), metaRGI$Timepoint, metaRGI$ID)
  colnames(myM)<-c("counts", "timePoint", "ID")
  Model<-lme(counts ~ poly(timePoint, 2), random = ~1 | ID, data = myM)
  rgipval[i]<-anova(Model)[2,4]
}
adjrgipval<-p.adjust(rgipval, method = "BH")
rgi<-c(sum(rgipval<0.05), sum(adjrgipval<0.05))

vsearchpval<-vector()
for (i in 1:nrow(vsearchT)) {
  myM<-data.frame(unlist(vsearchT[i, ]), metaVSEARCH$Timepoint, metaVSEARCH$ID)
  colnames(myM)<-c("counts", "timePoint", "ID")
  Model<-lme(counts ~ poly(timePoint, 2), random = ~1 | ID, data = myM)
  vsearchpval[i]<-anova(Model)[2,4]
}
adjvsearchpval<-p.adjust(vsearchpval, method = "BH")
vsearch<-c(sum(vsearchpval<0.05), sum(adjvsearchpval<0.05))

bigT<-rbind(amr, rgi, vsearch)
###########################

#AMRGeneFamily
#AMR
#keep only the rows present in amr counts table in master list
amrMatchT<-geneMatchT %>% filter(row_number() %in% match(rownames(amrT), geneMatchT$AMR.name))
amrMatchT<-na.omit(amrMatchT)
amrMatchT<-amrMatchT[, c(18,11)] #keep only the gene name and AMR classification columns

amrTypes<-unique(amrMatchT$AMR.Gene.Family) #all drug types
newamrT<-c()
for (i in 1:length(amrTypes)) {
  genes<-amrMatchT %>% filter(amrMatchT$AMR.Gene.Family == amrTypes[i])
  currentFrame<-amrT[genes$AMR.name, ]
  sumRow<-colSums(currentFrame)
  newamrT<-rbind(newamrT, sumRow)
}
rownames(newamrT)<-amrTypes
newamrT<-newamrT[, -which(colSums(newamrT) == 0)]

#RGI
#keep only the rows present in rgi counts table in master list
rgiMatchT<-geneMatchT %>% filter(row_number() %in% match(rownames(rgiT), geneMatchT$RGI.CARD.Short.Name))
rgiMatchT<-rgiMatchT[, c(15,11)] #keep only the gene name and AMR classification columns

rgiTypes<-unique(rgiMatchT$AMR.Gene.Family) #all drug types
newrgiT<-c()
for (i in 1:length(rgiTypes)) {
  genes<-rgiMatchT %>% filter(rgiMatchT$AMR.Gene.Family == rgiTypes[i])
  currentFrame<-rgiT[genes$RGI.CARD.Short.Name, ]
  sumRow<-colSums(currentFrame)
  newrgiT<-rbind(newrgiT, sumRow)
}
rownames(newrgiT)<-rgiTypes

#VSEARCH
#keep only the rows present in vsearch counts table in master list
vsearchMatchT<-geneMatchT %>% filter(row_number() %in% match(rownames(vsearchT), geneMatchT$Vsearch.ARO.Name))
vsearchMatchT<-vsearchMatchT[, c(7,11)] #keep only the gene name and AMR classification columns

vsearchTypes<-unique(vsearchMatchT$AMR.Gene.Family) #all drug types
newvsearchT<-c()
for (i in 1:length(vsearchTypes)) {
  genes<-vsearchMatchT %>% filter(vsearchMatchT$AMR.Gene.Family == vsearchTypes[i])
  currentFrame<-vsearchT[genes$Vsearch.ARO.Name, ]
  sumRow<-colSums(currentFrame)
  newvsearchT<-rbind(newvsearchT, sumRow)
}
rownames(newvsearchT)<-vsearchTypes

#meta for each table
metaAMR<-metaData[colnames(newamrT), ]
metaRGI<-metaData[colnames(newrgiT), ]
metaVSEARCH<-metaData[colnames(newvsearchT), ]

#amr 2nd order mixed lm anova pvals
AMRGeneFamily<-vector()
for (i in 1:nrow(newamrT)) {
  myM<-data.frame(unlist(newamrT[i, ]), metaAMR$Timepoint, metaAMR$ID)
  colnames(myM)<-c("counts", "timePoint", "ID")
  polyM<-lme(counts ~ poly(timePoint, 2), random = ~1 | ID, data = myM)
  
  AMRGeneFamily[i]<-anova(polyM)[2,4]
}
adjAMRGeneFamily<-p.adjust(AMRGeneFamily, method = "BH")
amr<-c(sum(AMRGeneFamily<0.05), sum(adjAMRGeneFamily<0.05))
#rgi 2nd order mixed lm anova pvals
newrgiT<-newrgiT[apply(newrgiT == 0, 1, sum) <= (ncol(newrgiT)*0.8), ]
AMRGeneFamily<-vector()
for (i in 1:nrow(newrgiT)) {
  myM<-data.frame(unlist(newrgiT[i, ]), metaRGI$Timepoint, metaRGI$ID)
  colnames(myM)<-c("counts", "timePoint", "ID")
  polyM<-lme(counts ~ poly(timePoint, 2), random = ~1 | ID, data = myM)
  
  AMRGeneFamily[i]<-anova(polyM)[2,4]
}
adjAMRGeneFamily<-p.adjust(AMRGeneFamily, method = "BH")
rgi<-c(sum(AMRGeneFamily<0.05), sum(adjAMRGeneFamily<0.05))
#vsearch 2nd order mixed lm anova pvals
AMRGeneFamily<-vector()
for (i in 1:nrow(newvsearchT)) {
  myM<-data.frame(unlist(newvsearchT[i, ]), metaVSEARCH$Timepoint, metaVSEARCH$ID)
  colnames(myM)<-c("counts", "timePoint", "ID")
  polyM<-lme(counts ~ poly(timePoint, 2), random = ~1 | ID, data = myM)
  
  AMRGeneFamily[i]<-anova(polyM)[2,4]
}
adjAMRGeneFamily<-p.adjust(AMRGeneFamily, method = "BH")
vsearch<-c(sum(AMRGeneFamily<0.05), sum(adjAMRGeneFamily<0.05))

bigT<-cbind(bigT, rbind(amr, rgi, vsearch))

#CorrectedAMRClassification
#AMR 
#keep only the rows present in amr counts table in master list
amrMatchT<-geneMatchT %>% filter(row_number() %in% match(rownames(amrT), geneMatchT$AMR.name))
amrMatchT<-na.omit(amrMatchT)
amrMatchT<-amrMatchT[, c(18,17)] #keep only the gene name and AMR classification columns

amrTypes<-unique(amrMatchT$Corrected.AMR.classification) #all drug types
newamrT<-c()
for (i in 1:length(amrTypes)) {
  genes<-amrMatchT %>% filter(amrMatchT$Corrected.AMR.classification == amrTypes[i])
  currentFrame<-amrT[genes$AMR.name, ]
  sumRow<-colSums(currentFrame)
  newamrT<-rbind(newamrT, sumRow)
}
rownames(newamrT)<-amrTypes
newamrT<-newamrT[, -which(colSums(newamrT) == 0)]

#RGI
#keep only the rows present in rgi counts table in master list
rgiMatchT<-geneMatchT %>% filter(row_number() %in% match(rownames(rgiT), geneMatchT$RGI.CARD.Short.Name))
rgiMatchT<-rgiMatchT[, c(15,17)] #keep only the gene name and AMR classification columns

rgiTypes<-unique(rgiMatchT$Corrected.AMR.classification) #all drug types
newrgiT<-c()
for (i in 1:length(rgiTypes)) {
  genes<-rgiMatchT %>% filter(rgiMatchT$Corrected.AMR.classification == rgiTypes[i])
  currentFrame<-rgiT[genes$RGI.CARD.Short.Name, ]
  sumRow<-colSums(currentFrame)
  newrgiT<-rbind(newrgiT, sumRow)
}
rownames(newrgiT)<-rgiTypes

#VSEARCH
#keep only the rows present in vsearch counts table in master list
vsearchMatchT<-geneMatchT %>% filter(row_number() %in% match(rownames(vsearchT), geneMatchT$Vsearch.ARO.Name))
vsearchMatchT<-vsearchMatchT[, c(7,17)] #keep only the gene name and AMR classification columns

vsearchTypes<-unique(vsearchMatchT$Corrected.AMR.classification) #all drug types
newvsearchT<-c()
for (i in 1:length(vsearchTypes)) {
  genes<-vsearchMatchT %>% filter(vsearchMatchT$Corrected.AMR.classification == vsearchTypes[i])
  currentFrame<-vsearchT[genes$Vsearch.ARO.Name, ]
  sumRow<-colSums(currentFrame)
  newvsearchT<-rbind(newvsearchT, sumRow)
}
rownames(newvsearchT)<-vsearchTypes

#meta for each table
metaAMR<-metaData[colnames(newamrT), ]
metaRGI<-metaData[colnames(newrgiT), ]
metaVSEARCH<-metaData[colnames(newvsearchT), ]

#amr 2nd order mixed lm anova pvals
CorrectedAMRClassification<-vector()
for (i in 1:nrow(newamrT)) {
  myM<-data.frame(unlist(newamrT[i, ]), metaAMR$Timepoint, metaAMR$ID)
  colnames(myM)<-c("counts", "timePoint", "ID")
  polyM<-lme(counts ~ poly(timePoint, 2), random = ~1 | ID, data = myM)
  
  CorrectedAMRClassification[i]<-anova(polyM)[2,4]
}
adjCorrectedAMRClassification<-p.adjust(CorrectedAMRClassification, method = "BH")
amr<-c(sum(CorrectedAMRClassification<0.05), sum(adjCorrectedAMRClassification<0.05))
#rgi 2nd order mixed lm anova pvals
CorrectedAMRClassification<-vector()
for (i in 1:nrow(newrgiT)) {
  myM<-data.frame(unlist(newrgiT[i, ]), metaRGI$Timepoint, metaRGI$ID)
  colnames(myM)<-c("counts", "timePoint", "ID")
  polyM<-lme(counts ~ poly(timePoint, 2), random = ~1 | ID, data = myM)
  
  CorrectedAMRClassification[i]<-anova(polyM)[2,4]
}
adjCorrectedAMRClassification<-p.adjust(CorrectedAMRClassification, method = "BH")
rgi<-c(sum(CorrectedAMRClassification<0.05), sum(adjCorrectedAMRClassification<0.05))
#vsearch 2nd order mixed lm anova pvals
CorrectedAMRClassification<-vector()
for (i in 1:nrow(newvsearchT)) {
  myM<-data.frame(unlist(newvsearchT[i, ]), metaVSEARCH$Timepoint, metaVSEARCH$ID)
  colnames(myM)<-c("counts", "timePoint", "ID")
  polyM<-lme(counts ~ poly(timePoint, 2), random = ~1 | ID, data = myM)
  
  CorrectedAMRClassification[i]<-anova(polyM)[2,4]
}
adjCorrectedAMRClassification<-p.adjust(CorrectedAMRClassification, method = "BH")
vsearch<-c(sum(CorrectedAMRClassification<0.05), sum(adjCorrectedAMRClassification<0.05))

bigT<-cbind(bigT, rbind(amr, rgi, vsearch))

#DrugClass
#AMR
#keep only the rows present in amr counts table in master list
amrMatchT<-geneMatchT %>% filter(row_number() %in% match(rownames(amrT), geneMatchT$AMR.name))
amrMatchT<-na.omit(amrMatchT)
amrMatchT<-amrMatchT[, c(18,12)] #keep only the gene name and AMR classification columns

amrTypes<-unique(amrMatchT$Drug.Class) #all drug types
newamrT<-c()
for (i in 1:length(amrTypes)) {
  genes<-amrMatchT %>% filter(amrMatchT$Drug.Class == amrTypes[i])
  currentFrame<-amrT[genes$AMR.name, ]
  sumRow<-colSums(currentFrame)
  newamrT<-rbind(newamrT, sumRow)
}
rownames(newamrT)<-amrTypes
newamrT<-newamrT[, -which(colSums(newamrT) == 0)]

#RGI
#keep only the rows present in rgi counts table in master list
rgiMatchT<-geneMatchT %>% filter(row_number() %in% match(rownames(rgiT), geneMatchT$RGI.CARD.Short.Name))
rgiMatchT<-rgiMatchT[, c(15,12)] #keep only the gene name and AMR classification columns

rgiTypes<-unique(rgiMatchT$Drug.Class) #all drug types
newrgiT<-c()
for (i in 1:length(rgiTypes)) {
  genes<-rgiMatchT %>% filter(rgiMatchT$Drug.Class == rgiTypes[i])
  currentFrame<-rgiT[genes$RGI.CARD.Short.Name, ]
  sumRow<-colSums(currentFrame)
  newrgiT<-rbind(newrgiT, sumRow)
}
rownames(newrgiT)<-rgiTypes

#VSEARCH
#keep only the rows present in vsearch counts table in master list
vsearchMatchT<-geneMatchT %>% filter(row_number() %in% match(rownames(vsearchT), geneMatchT$Vsearch.ARO.Name))
vsearchMatchT<-vsearchMatchT[, c(7,12)] #keep only the gene name and AMR classification columns

vsearchTypes<-unique(vsearchMatchT$Drug.Class) #all drug types
newvsearchT<-c()
for (i in 1:length(vsearchTypes)) {
  genes<-vsearchMatchT %>% filter(vsearchMatchT$Drug.Class == vsearchTypes[i])
  currentFrame<-vsearchT[genes$Vsearch.ARO.Name, ]
  sumRow<-colSums(currentFrame)
  newvsearchT<-rbind(newvsearchT, sumRow)
}
rownames(newvsearchT)<-vsearchTypes

#meta for each table
metaAMR<-metaData[colnames(newamrT), ]
metaRGI<-metaData[colnames(newrgiT), ]
metaVSEARCH<-metaData[colnames(newvsearchT), ]

#amr 2nd order mixed lm anova pvals
DrugClass<-vector()
for (i in 1:nrow(newamrT)) {
  myM<-data.frame(unlist(newamrT[i, ]), metaAMR$Timepoint, metaAMR$ID)
  colnames(myM)<-c("counts", "timePoint", "ID")
  polyM<-lme(counts ~ poly(timePoint, 2), random = ~1 | ID, data = myM)
  
  DrugClass[i]<-anova(polyM)[2,4]
}
adjDrugClass<-p.adjust(DrugClass, method = "BH")
amr<-c(sum(DrugClass<0.05), sum(adjDrugClass<0.05))
#rgi 2nd order mixed lm anova pvals
DrugClass<-vector()
for (i in 1:nrow(newrgiT)) {
  myM<-data.frame(unlist(newrgiT[i, ]), metaRGI$Timepoint, metaRGI$ID)
  colnames(myM)<-c("counts", "timePoint", "ID")
  polyM<-lme(counts ~ poly(timePoint, 2), random = ~1 | ID, data = myM)
  
  DrugClass[i]<-anova(polyM)[2,4]
}
adjDrugClass<-p.adjust(DrugClass, method = "BH")
rgi<-c(sum(DrugClass<0.05), sum(adjDrugClass<0.05))
#vsearch 2nd order mixed lm anova pvals
DrugClass<-vector()
for (i in 1:nrow(newvsearchT)) {
  myM<-data.frame(unlist(newvsearchT[i, ]), metaVSEARCH$Timepoint, metaVSEARCH$ID)
  colnames(myM)<-c("counts", "timePoint", "ID")
  polyM<-lme(counts ~ poly(timePoint, 2), random = ~1 | ID, data = myM)
  
  DrugClass[i]<-anova(polyM)[2,4]
}
adjDrugClass<-p.adjust(DrugClass, method = "BH")
vsearch<-c(sum(DrugClass<0.05), sum(adjDrugClass<0.05))

bigT<-cbind(bigT, rbind(amr, rgi, vsearch))

colnames(bigT)<-c("GenesOnly", "adjGenesOnly", "AMRGeneFamily", "adjAMRGeneFamily", "CorrectedAMRClassification", "adjCorrectedAMRClassification", "DrugClass", "adjDrugClass")

write.csv(bigT, "HitsTable.csv")
