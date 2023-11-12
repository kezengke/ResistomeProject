#filtering counts tables
rm(list = ls())

#counts tables
brackenT<-read.csv("CountsTables/brackenNormalized.csv", header = T, row.names = 1, check.names = F)
amrT<-read.csv("CountsTables/amrNormalized.csv", header = T, row.names = 1, check.names = F)
rgiT<-read.csv("CountsTables/rgiNormalized.csv", header = T, row.names = 1, check.names = F)
vsearchT<-read.csv("CountsTables/vsearchNormalized.csv", header = T, row.names = 1, check.names = F)
genusT<-read.csv("CountsTables/genusNormalized.csv", header = T, row.names = 1, check.names = F)
phylumT<-read.csv("CountsTables/phylumNormalized.csv", header = T, row.names = 1, check.names = F)
pathT<-read.csv("CountsTables/pathNormalized.csv", header = T, row.names = 1, check.names = F)

lowAbundance<-which(rowMeans(brackenT)<2)
brackenT<-brackenT[-lowAbundance, , drop = F]

lowAbundance<-which(rowMeans(genusT)<2)
genusT<-genusT[-lowAbundance, , drop = F]

lowAbundance<-which(rowMeans(phylumT)<2)
phylumT<-phylumT[-lowAbundance, , drop = F]

amrT<-amrT[apply(amrT == 0, 1, sum) <= (ncol(amrT)*0.8), ]
rgiT<-rgiT[apply(rgiT == 0, 1, sum) <= (ncol(rgiT)*0.8), ]
vsearchT<-vsearchT[apply(vsearchT == 0, 1, sum) <= (ncol(vsearchT)*0.8), ]

lowAbundance<-which(rowMeans(pathT)<2)
pathT<-pathT[-lowAbundance, , drop = F]

write.csv(brackenT, "CountsTables/brackenFiltered.csv")
write.csv(amrT, "CountsTables/amrFiltered.csv")
write.csv(rgiT, "CountsTables/rgiFiltered.csv")
write.csv(vsearchT, "CountsTables/vsearchFiltered.csv")
write.csv(genusT, "CountsTables/genusFiltered.csv")
write.csv(phylumT, "CountsTables/phylumFiltered.csv")
write.csv(pathT, "CountsTables/pathFiltered.csv")
