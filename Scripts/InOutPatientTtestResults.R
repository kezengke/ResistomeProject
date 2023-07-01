#t-test for in/out patients
rm(list = ls())

metaData<-read.csv("metaWithBins.csv", header = T, row.names = 1)

#counts tables
brackenT<-read.csv("CountsTables/brackenFiltered.csv", header = T, row.names = 1, check.names = F)
amrT<-read.csv("CountsTables/amrFiltered.csv", header = T, row.names = 1, check.names = F)
rgiT<-read.csv("CountsTables/rgiFiltered.csv", header = T, row.names = 1, check.names = F)
vsearchT<-read.csv("CountsTables/vsearchFiltered.csv", header = T, row.names = 1, check.names = F)
genusT<-read.csv("CountsTables/genusFiltered.csv", header = T, row.names = 1, check.names = F)

patients<-unique(metaData$ID)

pres<-vector()
posts<-vector()
for (i in 1:length(patients)) {
  patient<-metaData[metaData$ID == patients[i], , drop = F]
  pres[i]<-rownames(patient)[patient$bins == "PRE"]
  posts[i]<-rownames(patient)[which.min(abs(patient$Timepoint - 28))]
}

brackenPre<-brackenT[, pres, drop = F]
brackenPost<-brackenT[, posts, drop = F]
amrPre<-amrT[, pres, drop = F]
amrPost<-amrT[, intersect(colnames(amrT), posts), drop = F]
rgiPre<-rgiT[, pres, drop = F]
rgiPost<-rgiT[, posts, drop = F]
vsearchPre<-vsearchT[, pres, drop = F]
vsearchPost<-vsearchT[, posts, drop = F]
genusPre<-genusT[, pres, drop = F]
genusPost<-genusT[, posts, drop = F]

metaBRACKENpre<-metaData[colnames(brackenPre), , drop = F]
metaBRACKENpost<-metaData[colnames(brackenPost), , drop = F]
metaAMRpre<-metaData[colnames(amrPre), , drop = F]
metaAMRpost<-metaData[colnames(amrPost), , drop = F]
metaRGIpre<-metaData[colnames(rgiPre), , drop = F]
metaRGIpost<-metaData[colnames(rgiPost), , drop = F]
metaVSEARCHpre<-metaData[colnames(vsearchPre), , drop = F]
metaVSEARCHpost<-metaData[colnames(vsearchPost), , drop = F]
metaGENUSpre<-metaData[colnames(genusPre), , drop = F]
metaGENUSpost<-metaData[colnames(genusPost), , drop = F]

#p-values
brackenPreTtest_p<-apply(brackenPre, 1, function(x){t.test(unlist(x)~metaBRACKENpre$ptInOut)$p.value})
brackenPostTtest_p<-apply(brackenPost, 1, function(x){t.test(unlist(x)~metaBRACKENpost$ptInOut)$p.value})
amrPreTtest_p<-apply(amrPre, 1, function(x){t.test(unlist(x)~metaAMRpre$ptInOut)$p.value})
amrPostTtest_p<-apply(amrPost, 1, function(x){t.test(unlist(x)~metaAMRpost$ptInOut)$p.value})
rgiPreTtest_p<-apply(rgiPre, 1, function(x){t.test(unlist(x)~metaRGIpre$ptInOut)$p.value})
rgiPostTtest_p<-apply(rgiPost, 1, function(x){t.test(unlist(x)~metaRGIpost$ptInOut)$p.value})
vsearchPreTtest_p<-apply(vsearchPre, 1, function(x){t.test(unlist(x)~metaVSEARCHpre$ptInOut)$p.value})
vsearchPostTtest_p<-apply(vsearchPost, 1, function(x){t.test(unlist(x)~metaVSEARCHpost$ptInOut)$p.value})
genusPreTtest_p<-apply(genusPre, 1, function(x){t.test(unlist(x)~metaGENUSpre$ptInOut)$p.value})
genusPostTtest_p<-apply(genusPost, 1, function(x){t.test(unlist(x)~metaGENUSpost$ptInOut)$p.value})

#adj.p-values
brackenPreTtest_adj_p<-p.adjust(brackenPreTtest_p, method = "BH")
brackenPostTtest_adj_p<-p.adjust(brackenPostTtest_p, method = "BH")
amrPreTtest_adj_p<-p.adjust(amrPreTtest_p, method = "BH")
amrPostTtest_adj_p<-p.adjust(amrPostTtest_p, method = "BH")
rgiPreTtest_adj_p<-p.adjust(rgiPreTtest_p, method = "BH")
rgiPostTtest_adj_p<-p.adjust(rgiPostTtest_p, method = "BH")
vsearchPreTtest_adj_p<-p.adjust(vsearchPreTtest_p, method = "BH")
vsearchPostTtest_adj_p<-p.adjust(vsearchPostTtest_p, method = "BH")
genusPreTtest_adj_p<-p.adjust(genusPreTtest_p, method = "BH")
genusPostTtest_adj_p<-p.adjust(genusPostTtest_p, method = "BH")


brackenPreResults<-data.frame(brackenPreTtest_p, brackenPreTtest_adj_p)
brackenPostResults<-data.frame(brackenPostTtest_p, brackenPostTtest_adj_p)
amrPreResults<-data.frame(amrPreTtest_p, amrPreTtest_adj_p)
amrPostResults<-data.frame(amrPostTtest_p, amrPostTtest_adj_p)
rgiPreResults<-data.frame(rgiPreTtest_p, rgiPreTtest_adj_p)
rgiPostResults<-data.frame(rgiPostTtest_p, rgiPostTtest_adj_p)
vsearchPreResults<-data.frame(vsearchPreTtest_p, vsearchPreTtest_adj_p)
vsearchPostResults<-data.frame(vsearchPostTtest_p, vsearchPostTtest_adj_p)
genusPreResults<-data.frame(genusPreTtest_p, genusPreTtest_adj_p)
genusPostResults<-data.frame(genusPostTtest_p, genusPostTtest_adj_p)


rownames(brackenPreResults)<-rownames(brackenT)
rownames(brackenPostResults)<-rownames(brackenT)
rownames(amrPreResults)<-rownames(amrT)
rownames(amrPostResults)<-rownames(amrT)
rownames(rgiPreResults)<-rownames(rgiT)
rownames(rgiPostResults)<-rownames(rgiT)
rownames(vsearchPreResults)<-rownames(vsearchT)
rownames(vsearchPostResults)<-rownames(vsearchT)
rownames(genusPreResults)<-rownames(genusT)
rownames(genusPostResults)<-rownames(genusT)

write.csv(brackenPreResults, "Ttest/brackenPreInOut.csv")
write.csv(brackenPostResults, "Ttest/brackenPostInOut.csv")
write.csv(amrPreResults, "Ttest/amrPreInOut.csv")
write.csv(amrPostResults, "Ttest/amrPostInOut.csv")
write.csv(rgiPreResults, "Ttest/rgiPreInOut.csv")
write.csv(rgiPostResults, "Ttest/rgiPostInOut.csv")
write.csv(vsearchPreResults, "Ttest/vsearchPreInOut.csv")
write.csv(vsearchPostResults, "Ttest/vsearchPostInOut.csv")
write.csv(genusPreResults, "Ttest/genusPreInOut.csv")
write.csv(genusPostResults, "Ttest/genusPostInOut.csv")



