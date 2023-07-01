#bray-curtis distance matrix
rm(list = ls())
library("vegan")

#counts tables
brackenT<-read.csv("CountsTables/brackenFiltered.csv", header = T, row.names = 1, check.names = F)
amrT<-read.csv("CountsTables/amrFiltered.csv", header = T, row.names = 1, check.names = F)
rgiT<-read.csv("CountsTables/rgiFiltered.csv", header = T, row.names = 1, check.names = F)
vsearchT<-read.csv("CountsTables/vsearchFiltered.csv", header = T, row.names = 1, check.names = F)

distanceMatrix<-function(countsT, dissimilarityMetric){
  distanceT = as.matrix(vegdist(t(countsT), method = dissimilarityMetric))
  # diag(myDissimilarityMatrix) = NA
  return(distanceT)
}

brackenM<-distanceMatrix(brackenT, "bray")
amrM<-distanceMatrix(amrT, "bray")
rgiM<-distanceMatrix(rgiT, "bray")
vsearchM<-distanceMatrix(vsearchT, "bray")

write.csv(brackenM, "distanceMatrix/brackenDistanceMatrix.csv")
write.csv(amrM, "distanceMatrix/amrDistanceMatrix.csv")
write.csv(rgiM, "distanceMatrix/rgiDistanceMatrix.csv")
write.csv(vsearchM, "distanceMatrix/vsearchDistanceMatrix.csv")