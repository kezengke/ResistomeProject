rm(list = ls())
library(stringr)
library("RColorBrewer")
library("vegan")
library("dplyr")
library("nlme")

#gene category match table
geneMatchT<-read.csv("MASTER_AMRlist_2023_03.csv", sep = ",", header = T)
#metadata
dukeSamples<-read.csv("Duke_samples_meta.csv", header = T, sep = ",")
rownames(dukeSamples)<-dukeSamples$Seq_sample
rownames(dukeSamples)<-sapply(str_split(rownames(dukeSamples), "_", n = 2), `[`, 2)
dukeSamples$ID<-paste0("D", gsub("pre","" , sapply(str_split(dukeSamples$sample, "D", n = 3), `[`, 2), ignore.case = T)) #adding sample ID
dukeSamples$ID<-gsub("npe", "", dukeSamples$ID, ignore.case = T)
dukeSamples<-dukeSamples %>% mutate(bins = case_when(between(dukeSamples$Timepoint, -100, -2) ~ "PRE",
                                                     between(dukeSamples$Timepoint, -3, 10) ~ "D7",
                                                     between(dukeSamples$Timepoint, 11, 17) ~ "D14",
                                                     between(dukeSamples$Timepoint, 18, 24) ~ "D21",
                                                     between(dukeSamples$Timepoint, 25, 45) ~ "D35",
                                                     between(dukeSamples$Timepoint, 46, 75) ~ "D60",
                                                     between(dukeSamples$Timepoint, 76, 965) ~ "D100"))

#read in AMR
amrT<-read.delim("CountsTables/AMR_counts.tsv", sep = "\t", header = T, row.names = 2)
amrT<-amrT[,-1]
amrT<-amrT[, grepl("^D", colnames(amrT))]
colnames(amrT)<-gsub(".amrfinder.txt", "", colnames(amrT)) #get rid of useless info
SampleWithDots<-sapply(str_count(colnames(amrT), "\\."), `[`, 1) == 2
colnames(amrT)[SampleWithDots]<-sub("\\.", "", colnames(amrT)[SampleWithDots]) #get rid of random "."s in sample names
colnames(amrT)<-sub("\\.", "-", colnames(amrT)) #replace . with - to match with metadata
colnames(amrT)<-sapply(str_split(colnames(amrT), "_", n = 2), `[`, 2)
amrT<-amrT[, -27]

#read in RGI
rgiT<-read.delim("CountsTables/RGI_counts.tsv", sep = "\t", header = T, row.names = 2)
rgiT<-rgiT[, -1]
rgiT<-rgiT[, grepl("^D", colnames(rgiT))]
colnames(rgiT)<-gsub(".rgi.txt", "", colnames(rgiT)) #get rid of useless info
SampleWithDots<-sapply(str_count(colnames(rgiT), "\\."), `[`, 1) == 2
colnames(rgiT)[SampleWithDots]<-sub("\\.", "", colnames(rgiT)[SampleWithDots]) #get rid of random "."s in sample names
colnames(rgiT)<-sub("\\.", "-", colnames(rgiT)) #replace . with - to match with metadata
colnames(rgiT)<-sapply(str_split(colnames(rgiT), "_", n = 2), `[`, 2)

#read in vsearch
vsearchT<-read.delim("CountsTables/vsearch_counts.tsv", sep = "\t", header = T, row.names = 2)
vsearchT<-vsearchT[, -1] #get rid of index column
vsearchT<-vsearchT[, grepl("^D", colnames(vsearchT))]
colnames(vsearchT)<-gsub(".txt", "", colnames(vsearchT)) #get rid of useless info

SampleWithDots<-sapply(str_count(colnames(vsearchT), "\\."), `[`, 1) == 2
colnames(vsearchT)[SampleWithDots]<-sub("\\.", "", colnames(vsearchT)[SampleWithDots]) #get rid of random "."s in sample names
colnames(vsearchT)<-sub("\\.", "-", colnames(vsearchT)) #replace . with - to match with metadata
colnames(vsearchT)<-sapply(str_split(colnames(vsearchT), "_", n = 2), `[`, 2)

rownames(vsearchT)<-sapply(str_split(rownames(vsearchT), "\\|", n = 6), `[`, 6)


bigT<-c()
###########################
#Normalization
n<-colSums(amrT)
sumx<-sum(amrT)
normamrT<-amrT
for (i in 1:ncol(normamrT)) {
  normamrT[,i]<-normamrT[,i]/n[i]
}
normamrT<-log10(normamrT*(sumx/ncol(normamrT))+1)

n<-colSums(rgiT)
sumx<-sum(rgiT)
normrgiT<-rgiT
for (i in 1:ncol(normrgiT)) {
  normrgiT[,i]<-normrgiT[,i]/n[i]
}
normrgiT<-log10(normrgiT*(sumx/ncol(normrgiT))+1)

n<-colSums(vsearchT)
sumx<-sum(vsearchT)
normvsearchT<-vsearchT
for (i in 1:ncol(normvsearchT)) {
  normvsearchT[,i]<-normvsearchT[,i]/n[i]
}
normvsearchT<-log10(normvsearchT*(sumx/ncol(normvsearchT))+1)

#meta for each table
metaAMR<-dukeSamples[colnames(normamrT), ]
metaRGI<-dukeSamples[colnames(normrgiT), ]
metaVsearch<-dukeSamples[colnames(normvsearchT), ]

#lm 2nd order
normamrT<-normamrT[apply(normamrT == 0, 1, sum) <= (ncol(normamrT)*0.8), ]
amrpval<-vector()
for (i in 1:nrow(normamrT)) {
  myM<-data.frame(unlist(normamrT[i, ]), metaAMR$Timepoint, metaAMR$ID)
  colnames(myM)<-c("counts", "timePoint", "ID")
  Model<-lme(counts ~ poly(timePoint, 2), random = ~1 | ID, data = myM)
  amrpval[i]<-anova(Model)[2,4]
}
adjamrpval<-p.adjust(amrpval, method = "BH")
amr<-c(sum(amrpval<0.05), sum(adjamrpval<0.05))

normrgiT<-normrgiT[apply(normrgiT == 0, 1, sum) <= (ncol(normrgiT)*0.8), ]
rgipval<-vector()
for (i in 1:nrow(normrgiT)) {
  myM<-data.frame(unlist(normrgiT[i, ]), metaRGI$Timepoint, metaRGI$ID)
  colnames(myM)<-c("counts", "timePoint", "ID")
  Model<-lme(counts ~ poly(timePoint, 2), random = ~1 | ID, data = myM)
  rgipval[i]<-anova(Model)[2,4]
}
adjrgipval<-p.adjust(rgipval, method = "BH")
rgi<-c(sum(rgipval<0.05), sum(adjrgipval<0.05))

normvsearchT<-normvsearchT[apply(normvsearchT == 0, 1, sum) <= (ncol(normvsearchT)*0.8), ]
vsearchpval<-vector()
for (i in 1:nrow(normvsearchT)) {
  myM<-data.frame(unlist(normvsearchT[i, ]), metaVsearch$Timepoint, metaVsearch$ID)
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
metaAMR<-dukeSamples[colnames(newamrT), ]
metaRGI<-dukeSamples[colnames(newrgiT), ]
metaVsearch<-dukeSamples[colnames(newvsearchT), ]

#Normalization
n<-colSums(newamrT)
sumx<-sum(newamrT)
for (i in 1:ncol(newamrT)) {
  newamrT[,i]<-newamrT[,i]/n[i]
}
newamrT<-log10(newamrT*(sumx/ncol(newamrT))+1)

n<-colSums(newrgiT)
sumx<-sum(newrgiT)
for (i in 1:ncol(newrgiT)) {
  newrgiT[,i]<-newrgiT[,i]/n[i]
}
newrgiT<-log10(newrgiT*(sumx/ncol(newrgiT))+1)

n<-colSums(newvsearchT)
sumx<-sum(newvsearchT)
for (i in 1:ncol(newvsearchT)) {
  newvsearchT[,i]<-newvsearchT[,i]/n[i]
}
newvsearchT<-log10(newvsearchT*(sumx/ncol(newvsearchT))+1)

#amr 2nd order mixed lm anova pvals
newamrT<-newamrT[apply(newamrT == 0, 1, sum) <= (ncol(newamrT)*0.8), ]
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
newvsearchT<-newvsearchT[apply(newvsearchT == 0, 1, sum) <= (ncol(newvsearchT)*0.8), ]
AMRGeneFamily<-vector()
for (i in 1:nrow(newvsearchT)) {
  myM<-data.frame(unlist(newvsearchT[i, ]), metaVsearch$Timepoint, metaVsearch$ID)
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
metaAMR<-dukeSamples[colnames(newamrT), ]
metaRGI<-dukeSamples[colnames(newrgiT), ]
metaVsearch<-dukeSamples[colnames(newvsearchT), ]

#Normalization
n<-colSums(newamrT)
sumx<-sum(newamrT)
for (i in 1:ncol(newamrT)) {
  newamrT[,i]<-newamrT[,i]/n[i]
}
newamrT<-log10(newamrT*(sumx/ncol(newamrT))+1)

n<-colSums(newrgiT)
sumx<-sum(newrgiT)
for (i in 1:ncol(newrgiT)) {
  newrgiT[,i]<-newrgiT[,i]/n[i]
}
newrgiT<-log10(newrgiT*(sumx/ncol(newrgiT))+1)

n<-colSums(newvsearchT)
sumx<-sum(newvsearchT)
for (i in 1:ncol(newvsearchT)) {
  newvsearchT[,i]<-newvsearchT[,i]/n[i]
}
newvsearchT<-log10(newvsearchT*(sumx/ncol(newvsearchT))+1)

#amr 2nd order mixed lm anova pvals
newamrT<-newamrT[apply(newamrT == 0, 1, sum) <= (ncol(newamrT)*0.8), ]
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
newrgiT<-newrgiT[apply(newrgiT == 0, 1, sum) <= (ncol(newrgiT)*0.8), ]
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
newvsearchT<-newvsearchT[apply(newvsearchT == 0, 1, sum) <= (ncol(newvsearchT)*0.8), ]
CorrectedAMRClassification<-vector()
for (i in 1:nrow(newvsearchT)) {
  myM<-data.frame(unlist(newvsearchT[i, ]), metaVsearch$Timepoint, metaVsearch$ID)
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
metaAMR<-dukeSamples[colnames(newamrT), ]
metaRGI<-dukeSamples[colnames(newrgiT), ]
metaVsearch<-dukeSamples[colnames(newvsearchT), ]

#Normalization
n<-colSums(newamrT)
sumx<-sum(newamrT)
for (i in 1:ncol(newamrT)) {
  newamrT[,i]<-newamrT[,i]/n[i]
}
newamrT<-log10(newamrT*(sumx/ncol(newamrT))+1)

n<-colSums(newrgiT)
sumx<-sum(newrgiT)
for (i in 1:ncol(newrgiT)) {
  newrgiT[,i]<-newrgiT[,i]/n[i]
}
newrgiT<-log10(newrgiT*(sumx/ncol(newrgiT))+1)

n<-colSums(newvsearchT)
sumx<-sum(newvsearchT)
for (i in 1:ncol(newvsearchT)) {
  newvsearchT[,i]<-newvsearchT[,i]/n[i]
}
newvsearchT<-log10(newvsearchT*(sumx/ncol(newvsearchT))+1)

#amr 2nd order mixed lm anova pvals
newamrT<-newamrT[apply(newamrT == 0, 1, sum) <= (ncol(newamrT)*0.8), ]
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
newrgiT<-newrgiT[apply(newrgiT == 0, 1, sum) <= (ncol(newrgiT)*0.8), ]
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
newvsearchT<-newvsearchT[apply(newvsearchT == 0, 1, sum) <= (ncol(newvsearchT)*0.8), ]
DrugClass<-vector()
for (i in 1:nrow(newvsearchT)) {
  myM<-data.frame(unlist(newvsearchT[i, ]), metaVsearch$Timepoint, metaVsearch$ID)
  colnames(myM)<-c("counts", "timePoint", "ID")
  polyM<-lme(counts ~ poly(timePoint, 2), random = ~1 | ID, data = myM)
  
  DrugClass[i]<-anova(polyM)[2,4]
}
adjDrugClass<-p.adjust(DrugClass, method = "BH")
vsearch<-c(sum(DrugClass<0.05), sum(adjDrugClass<0.05))

bigT<-cbind(bigT, rbind(amr, rgi, vsearch))

colnames(bigT)<-c("GenesOnly", "adjGenesOnly", "AMRGeneFamily", "adjAMRGeneFamily", "CorrectedAMRClassification", "adjCorrectedAMRClassification", "DrugClass", "adjDrugClass")

write.table(bigT, "HitsTable.txt", sep = "\t")
