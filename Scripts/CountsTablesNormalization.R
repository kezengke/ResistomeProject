#normalizing counts tables
rm(list = ls())

#counts tables
brackenT<-read.csv("CountsTables/brackenProcessed.csv", header = T, row.names = 1, check.names = F)
amrT<-read.csv("CountsTables/amrProcessed.csv", header = T, row.names = 1, check.names = F)
rgiT<-read.csv("CountsTables/rgiProcessed.csv", header = T, row.names = 1, check.names = F)
vsearchT<-read.csv("CountsTables/vsearchProcessed.csv", header = T, row.names = 1, check.names = F)
genusT<-read.csv("CountsTables/genusProcessed.csv", header = T, row.names = 1, check.names = F)

#Normalization
n<-colSums(brackenT)
sumx<-sum(brackenT)
for (i in 1:ncol(brackenT)) {
  brackenT[,i]<-brackenT[,i]/n[i]
}
brackenT<-log10(brackenT*(sumx/ncol(brackenT))+1)

n<-colSums(amrT)
sumx<-sum(amrT)
for (i in 1:ncol(amrT)) {
  amrT[,i]<-amrT[,i]/n[i]
}
amrT<-log10(amrT*(sumx/ncol(amrT))+1)

n<-colSums(rgiT)
sumx<-sum(rgiT)
for (i in 1:ncol(rgiT)) {
  rgiT[,i]<-rgiT[,i]/n[i]
}
rgiT<-log10(rgiT*(sumx/ncol(rgiT))+1)

n<-colSums(vsearchT)
sumx<-sum(vsearchT)
for (i in 1:ncol(vsearchT)) {
  vsearchT[,i]<-vsearchT[,i]/n[i]
}
vsearchT<-log10(vsearchT*(sumx/ncol(vsearchT))+1)

n<-colSums(genusT)
sumx<-sum(genusT)
for (i in 1:ncol(genusT)) {
  genusT[,i]<-genusT[,i]/n[i]
}
genusT<-log10(genusT*(sumx/ncol(genusT))+1)

write.csv(brackenT, "CountsTables/brackenNormalized.csv")
write.csv(amrT, "CountsTables/amrNormalized.csv")
write.csv(rgiT, "CountsTables/rgiNormalized.csv")
write.csv(vsearchT, "CountsTables/vsearchNormalized.csv")
write.csv(genusT, "CountsTables/genusNormalized.csv")
