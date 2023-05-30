#wilcoxon test for in/out patients
rm(list = ls())

metaData<-read.csv("metaWithBins.csv", header = T, row.names = 1)

#counts tables
brackenT<-read.csv("CountsTables/brackenFiltered.csv", header = T, row.names = 1, check.names = F)
amrT<-read.csv("CountsTables/amrFiltered.csv", header = T, row.names = 1, check.names = F)
rgiT<-read.csv("CountsTables/rgiFiltered.csv", header = T, row.names = 1, check.names = F)
vsearchT<-read.csv("CountsTables/vsearchFiltered.csv", header = T, row.names = 1, check.names = F)

pres<-rownames(metaData)[metaData$bins == "PRE"]
posts<-rownames(metaData)[metaData$bins != "PRE"]

brackenPre<-brackenT[, pres, drop = F]
brackenPost<-brackenT[, posts, drop = F]
amrPre<-amrT[, pres, drop = F]
amrPost<-amrT[, intersect(colnames(amrT), posts), drop = F]
rgiPre<-rgiT[, pres, drop = F]
rgiPost<-rgiT[, posts, drop = F]
vsearchPre<-vsearchT[, pres, drop = F]
vsearchPost<-vsearchT[, posts, drop = F]

metaBRACKENpre<-metaData[colnames(brackenPre), , drop = F]
metaBRACKENpost<-metaData[colnames(brackenPost), , drop = F]
metaAMRpre<-metaData[colnames(amrPre), , drop = F]
metaAMRpost<-metaData[colnames(amrPost), , drop = F]
metaRGIpre<-metaData[colnames(rgiPre), , drop = F]
metaRGIpost<-metaData[colnames(rgiPost), , drop = F]
metaVSEARCHpre<-metaData[colnames(vsearchPre), , drop = F]
metaVSEARCHpost<-metaData[colnames(vsearchPost), , drop = F]

#p-values
brackenPreWilcox_p<-apply(brackenPre, 1, function(x){wilcox.test(unlist(x)~metaBRACKENpre$ptInOut)$p.value})
brackenPostWilcox_p<-apply(brackenPost, 1, function(x){wilcox.test(unlist(x)~metaBRACKENpost$ptInOut)$p.value})
amrPreWilcox_p<-apply(amrPre, 1, function(x){wilcox.test(unlist(x)~metaAMRpre$ptInOut)$p.value})
amrPostWilcox_p<-apply(amrPost, 1, function(x){wilcox.test(unlist(x)~metaAMRpost$ptInOut)$p.value})
rgiPreWilcox_p<-apply(rgiPre, 1, function(x){wilcox.test(unlist(x)~metaRGIpre$ptInOut)$p.value})
rgiPostWilcox_p<-apply(rgiPost, 1, function(x){wilcox.test(unlist(x)~metaRGIpost$ptInOut)$p.value})
vsearchPreWilcox_p<-apply(vsearchPre, 1, function(x){wilcox.test(unlist(x)~metaVSEARCHpre$ptInOut)$p.value})
vsearchPostWilcox_p<-apply(vsearchPost, 1, function(x){wilcox.test(unlist(x)~metaVSEARCHpost$ptInOut)$p.value})

#adj.p-values
brackenPreWilcox_adj_p<-p.adjust(brackenPreWilcox_p, method = "BH")
brackenPostWilcox_adj_p<-p.adjust(brackenPostWilcox_p, method = "BH")
amrPreWilcox_adj_p<-p.adjust(amrPreWilcox_p, method = "BH")
amrPostWilcox_adj_p<-p.adjust(amrPostWilcox_p, method = "BH")
rgiPreWilcox_adj_p<-p.adjust(rgiPreWilcox_p, method = "BH")
rgiPostWilcox_adj_p<-p.adjust(rgiPostWilcox_p, method = "BH")
vsearchPreWilcox_adj_p<-p.adjust(vsearchPreWilcox_p, method = "BH")
vsearchPostWilcox_adj_p<-p.adjust(vsearchPostWilcox_p, method = "BH")


brackenPreResults<-data.frame(brackenPreWilcox_p, brackenPreWilcox_adj_p)
brackenPostResults<-data.frame(brackenPostWilcox_p, brackenPostWilcox_adj_p)
amrPreResults<-data.frame(amrPreWilcox_p, amrPreWilcox_adj_p)
amrPostResults<-data.frame(amrPostWilcox_p, amrPostWilcox_adj_p)
rgiPreResults<-data.frame(rgiPreWilcox_p, rgiPreWilcox_adj_p)
rgiPostResults<-data.frame(rgiPostWilcox_p, rgiPostWilcox_adj_p)
vsearchPreResults<-data.frame(vsearchPreWilcox_p, vsearchPreWilcox_adj_p)
vsearchPostResults<-data.frame(vsearchPostWilcox_p, vsearchPostWilcox_adj_p)


rownames(brackenPreResults)<-rownames(brackenT)
rownames(brackenPostResults)<-rownames(brackenT)
rownames(amrPreResults)<-rownames(amrT)
rownames(amrPostResults)<-rownames(amrT)
rownames(rgiPreResults)<-rownames(rgiT)
rownames(rgiPostResults)<-rownames(rgiT)
rownames(vsearchPreResults)<-rownames(vsearchT)
rownames(vsearchPostResults)<-rownames(vsearchT)


write.csv(brackenPreResults, "Wilcox/brackenPreInOut.csv")
write.csv(brackenPostResults, "Wilcox/brackenPostInOut.csv")
write.csv(amrPreResults, "Wilcox/amrPreInOut.csv")
write.csv(amrPostResults, "Wilcox/amrPostInOut.csv")
write.csv(rgiPreResults, "Wilcox/rgiPreInOut.csv")
write.csv(rgiPostResults, "Wilcox/rgiPostInOut.csv")
write.csv(vsearchPreResults, "Wilcox/vsearchPreInOut.csv")
write.csv(vsearchPostResults, "Wilcox/vsearchPostInOut.csv")




