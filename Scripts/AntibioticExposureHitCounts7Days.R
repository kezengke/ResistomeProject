rm(list = ls())
library("dplyr")
library(stringr)

A_expo<-read.csv("AntibioticsExposure.csv", check.names = F, header = T)

#gene category match table
geneMatchT<-read.csv("MASTER_AMRlist_2023_03.csv", sep = ",", header = T)

#metadata
dukeSamples<-read.csv("Duke_samples_meta.csv", header = T, sep = ",")
rownames(dukeSamples)<-dukeSamples$Seq_sample
rownames(dukeSamples)<-sapply(str_split(rownames(dukeSamples), "_", n = 2), `[`, 2)
dukeSamples$ID<-paste0("D", gsub("pre","" , sapply(str_split(dukeSamples$sample, "D", n = 3), `[`, 2), ignore.case = T)) #adding sample ID
dukeSamples$ID<-gsub("npe", "", dukeSamples$ID, ignore.case = T)

IDtype<-unique(A_expo$ID)
ABtype<-unique(A_expo$AntibioticType)

metaWithAB<-dukeSamples
for (i in 1:length(ABtype)) {
  AB<-vector()
  checkABframe<-A_expo[A_expo$AntibioticType == ABtype[i], ]
  for (j in 1:nrow(dukeSamples)) {
    AB[j]<-sum(dukeSamples[j, 2] > checkABframe[checkABframe$ID == dukeSamples[j, 1], 3] & 
                 dukeSamples[j, 2] < (checkABframe[checkABframe$ID == dukeSamples[j, 1], 4]+7))
  }
  new_AB<-ifelse(AB > 0, "Yes", "No")
  cname<-ABtype[i]
  metaWithAB[, cname]<-new_AB
}

#bracken
brackenT<-read.delim("CountsTables/duke_bracken.csv", sep = ",", header = T, row.names = 1)
brackenT<-brackenT[, -c(1,2)] #get rid of taxonomy ID and level
brackenT<-brackenT[, grepl("num", colnames(brackenT))] #get rid of fractions

colnames(brackenT)<-gsub(".bracken.out_num", "", colnames(brackenT)) #get rid of useless info
SampleWithDots<-sapply(str_count(colnames(brackenT), "\\."), `[`, 1) == 2
colnames(brackenT)[SampleWithDots]<-sub("\\.", "", colnames(brackenT)[SampleWithDots]) #get rid of random "."s in sample names
colnames(brackenT)<-sub("\\.", "-", colnames(brackenT)) #replace . with - to match with metadata
colnames(brackenT)<-sapply(str_split(colnames(brackenT), "_", n = 2), `[`, 2) #keeping sequencing info for matching

#Normalization
n<-colSums(brackenT)
sumx<-sum(brackenT)
for (i in 1:ncol(brackenT)) {
  brackenT[,i]<-brackenT[,i]/n[i]
}
brackenT<-log10(brackenT*(sumx/ncol(brackenT))+1)
#filter
brackenT<-brackenT[apply(brackenT == 0, 1, sum) <= (ncol(brackenT)*0.8), ]

#meta for each table
metaBracken<-metaWithAB[colnames(brackenT), 6:ncol(metaWithAB)]
#exclude antibiotics that are not exposed to any samples

metaBracken<-metaBracken[, -which(colSums(metaBracken == "Yes") < 2)]

hits<-vector()
for (i in 1:ncol(metaBracken)) {
  testingMeta<-metaBracken[, i, drop = F]
  pval<-vector()
  for (j in 1:nrow(brackenT)) {
    pval[j]<-t.test(unlist(brackenT[j, ])~factor(testingMeta[,1]))$p.value
  }
  hits[i]<-sum(pval<0.05)
}
