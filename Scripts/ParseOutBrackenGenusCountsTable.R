#Parse out genus level from bracken species counts table
rm(list = ls())
library(stringr)

brackenT<-read.csv("CountsTables/brackenProcessed.csv", header = T, row.names = 1, check.names = F)
rownames(brackenT)<-gsub("\\[|\\]", "", rownames(brackenT))
genus<-rownames(brackenT)
genus<-sapply(str_split(genus, " ", n = 2), `[`, 1)

genusType<-unique(genus)

genusT<-vector()
for (taxa in genusType) {
  allRows<-which(genus == taxa)
  currentTaxa<-colSums(brackenT[allRows, , drop = F])
  genusT<-rbind(genusT, currentTaxa)
}
rownames(genusT)<-genusType

write.csv(genusT, "CountsTables/genusProcessed.csv")
