#t-test for in/out patients
rm(list = ls())

metaData<-read.csv("MetaBinnedWithInOutDonorTreatmentType.csv", header = T, row.names = 1)

#counts tables
speciesT<-read.csv("CountsTables/brackenFiltered.csv", header = T, row.names = 1, check.names = F)
genusT<-read.csv("CountsTables/genusFiltered.csv", header = T, row.names = 1, check.names = F)
phylumT<-read.csv("CountsTables/phylumFiltered.csv", header = T, row.names = 1, check.names = F)
amrT<-read.csv("CountsTables/amrFiltered.csv", header = T, row.names = 1, check.names = F)
rgiT<-read.csv("CountsTables/rgiFiltered.csv", header = T, row.names = 1, check.names = F)
vsearchT<-read.csv("CountsTables/vsearchFiltered.csv", header = T, row.names = 1, check.names = F)
pathwayT<-read.csv("CountsTables/pathFiltered.csv", header = T, row.names = 1, check.names = F)

patients<-unique(metaData$ID)

pres<-rownames(metaData)[metaData$bins == "PRE"]
posts<-rownames(metaData)[metaData$bins == "D30"]

ttestFun <- function(counts, pres, posts, meta, tableType) {
  counts<-counts[, intersect(colnames(counts), rownames(meta)), drop = F]
  meta<-meta[intersect(colnames(counts), rownames(meta)), , drop = F]
  
  metaPost<-meta[posts, , drop = F]
  
  metaPre<-meta[pres, , drop = F]
  metaPre<-metaPre[metaPre$ID %in% metaPost$ID, , drop = F]
  
  PreT<-counts[, rownames(metaPre), drop = F]
  PostT<-counts[, rownames(metaPost), drop = F]
  
  #p-values
  PreTtest_p<-apply(PreT, 1, function(x){t.test(unlist(x)~metaPre$ptInOut)$p.value})
  PostTtest_p<-apply(PostT, 1, function(x){t.test(unlist(x)~metaPost$ptInOut)$p.value})

  #adj.p-values
  PreTtest_adj_p<-p.adjust(PreTtest_p, method = "BH")
  PostTtest_adj_p<-p.adjust(PostTtest_p, method = "BH")

  PreResults<-data.frame(PreTtest_p, PreTtest_adj_p)
  PostResults<-data.frame(PostTtest_p, PostTtest_adj_p)

  rownames(PreResults)<-rownames(counts)
  rownames(PostResults)<-rownames(counts)

  write.csv(PreResults, paste0("Ttest/", tableType, "PreInOut.csv"))
  write.csv(PostResults, paste0("Ttest/", tableType, "PostInOut.csv"))

}

ttestFun(counts = speciesT, pres = pres, posts = posts, meta = metaData, tableType = "species")
ttestFun(counts = genusT, pres = pres, posts = posts, meta = metaData, tableType = "genus")
ttestFun(counts = phylumT, pres = pres, posts = posts, meta = metaData, tableType = "phylum")

ttestFun(counts = amrT, pres = pres, posts = posts, meta = metaData, tableType = "AMR")
ttestFun(counts = rgiT, pres = pres, posts = posts, meta = metaData, tableType = "RGI")
ttestFun(counts = vsearchT, pres = pres, posts = posts, meta = metaData, tableType = "vsearch")

ttestFun(counts = pathwayT, pres = pres, posts = posts, meta = metaData, tableType = "pathway")
