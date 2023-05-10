rm(list = ls())
library("dplyr")
library(stringr)

A_expo<-read.csv("AntibioticsExposure.csv", check.names = F, header = T)

#gene category match table
geneMatchT<-read.csv("MASTER_AMRlist_2023_03.csv", sep = ",", header = T)

#metadata
meta7days<-read.csv("ABexposure/ABmeta7days.csv", header = T, row.names = 1, check.names = F)
# meta10days<-read.csv("ABexposure/ABmeta10days.csv", header = T, row.names = 1, check.names = F)
# meta30days<-read.csv("ABexposure/ABmeta30days.csv", header = T, row.names = 1, check.names = F)
# meta60days<-read.csv("ABexposure/ABmeta60days.csv", header = T, row.names = 1, check.names = F)
# meta90days<-read.csv("ABexposure/ABmeta90days.csv", header = T, row.names = 1, check.names = F)

#bracken
brackenT<-read.csv("CountsTables/brackenNormalized.csv", header = T, row.names = 1, check.names = F)
#filter
brackenT<-brackenT[apply(brackenT == 0, 1, sum) <= (ncol(brackenT)*0.8), ]

#meta for each table
metaBracken7<-meta7days[colnames(brackenT), 7:ncol(meta7days)]
#exclude antibiotics that are not exposed to any samples
# metaBracken7<-metaBracken7[, -which(colSums(metaBracken7 == "Yes") < 2)]
# metaBracken7<-metaBracken7[, -which(colSums(metaBracken7 == "No") < 2)]

hits7<-vector()
for (i in 1:ncol(metaBracken7)) {
  testingMeta<-metaBracken7[, i, drop = F]
  pval<-vector()
  for (j in 1:nrow(brackenT)) {
    if (colSums(testingMeta == "Yes") < 2) {
      pval[j]<-NA
    }
    else
      pval[j]<-t.test(unlist(brackenT[j, ])~factor(testingMeta[,1]))$p.value
  }
  hits7[i]<-sum(pval<0.05)
}
