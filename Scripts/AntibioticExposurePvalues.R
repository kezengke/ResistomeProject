rm(list = ls())
library("dplyr")
library(stringr)

A_expo<-read.csv("AntibioticsExposure.csv", check.names = F, header = T)

#gene category match table
geneMatchT<-read.csv("MASTER_AMRlist_2023_03.csv", sep = ",", header = T)

#metadata
meta7days<-read.csv("ABexposure/ABmeta7days.csv", header = T, row.names = 1, check.names = F)
meta10days<-read.csv("ABexposure/ABmeta10days.csv", header = T, row.names = 1, check.names = F)
meta30days<-read.csv("ABexposure/ABmeta30days.csv", header = T, row.names = 1, check.names = F)
meta60days<-read.csv("ABexposure/ABmeta60days.csv", header = T, row.names = 1, check.names = F)
meta90days<-read.csv("ABexposure/ABmeta90days.csv", header = T, row.names = 1, check.names = F)

#counts tables
brackenT<-read.csv("CountsTables/brackenFiltered.csv", header = T, row.names = 1, check.names = F)
amrT<-read.csv("CountsTables/amrFiltered.csv", header = T, row.names = 1, check.names = F)
rgiT<-read.csv("CountsTables/rgiFiltered.csv", header = T, row.names = 1, check.names = F)
vsearchT<-read.csv("CountsTables/vsearchFiltered.csv", header = T, row.names = 1, check.names = F)

#meta for each table
metaBRACKEN7<-meta7days[colnames(brackenT), 8:ncol(meta7days)]
metaBRACKEN10<-meta10days[colnames(brackenT), 8:ncol(meta10days)]
metaBRACKEN30<-meta30days[colnames(brackenT), 8:ncol(meta30days)]
metaBRACKEN60<-meta60days[colnames(brackenT), 8:ncol(meta60days)]
metaBRACKEN90<-meta90days[colnames(brackenT), 8:ncol(meta90days)]

metaAMR7<-meta7days[colnames(amrT), 8:ncol(meta7days)]
metaAMR10<-meta10days[colnames(amrT), 8:ncol(meta10days)]
metaAMR30<-meta30days[colnames(amrT), 8:ncol(meta30days)]
metaAMR60<-meta60days[colnames(amrT), 8:ncol(meta60days)]
metaAMR90<-meta90days[colnames(amrT), 8:ncol(meta90days)]

metaRGI7<-meta7days[colnames(rgiT), 8:ncol(meta7days)]
metaRGI10<-meta10days[colnames(rgiT), 8:ncol(meta10days)]
metaRGI30<-meta30days[colnames(rgiT), 8:ncol(meta30days)]
metaRGI60<-meta60days[colnames(rgiT), 8:ncol(meta60days)]
metaRGI90<-meta90days[colnames(rgiT), 8:ncol(meta90days)]

metaVSEARCH7<-meta7days[colnames(vsearchT), 8:ncol(meta7days)]
metaVSEARCH10<-meta10days[colnames(vsearchT), 8:ncol(meta10days)]
metaVSEARCH30<-meta30days[colnames(vsearchT), 8:ncol(meta30days)]
metaVSEARCH60<-meta60days[colnames(vsearchT), 8:ncol(meta60days)]
metaVSEARCH90<-meta90days[colnames(vsearchT), 8:ncol(meta90days)]

pvalCalc<-function(meta, countsT) {
  typepval<-sapply(meta, function(testingMeta) {
    if (sum(testingMeta == "Yes") < 2) {
      rep(NA, nrow(countsT))
    } else {
      apply(countsT, 1, function(row) {
        t.test(row~factor(testingMeta))$p.value
      })
    }
  })
  rownames(typepval) <- rownames(countsT)
  colnames(typepval) <- colnames(meta)
  return(typepval)
}

#calculating pvalues
bracken7<-pvalCalc(metaBRACKEN7, brackenT)
bracken10<-pvalCalc(metaBRACKEN10, brackenT)
bracken30<-pvalCalc(metaBRACKEN30, brackenT)
bracken60<-pvalCalc(metaBRACKEN60, brackenT)
bracken90<-pvalCalc(metaBRACKEN90, brackenT)

amr7<-pvalCalc(metaAMR7, amrT)
amr10<-pvalCalc(metaAMR10, amrT)
amr30<-pvalCalc(metaAMR30, amrT)
amr60<-pvalCalc(metaAMR60, amrT)
amr90<-pvalCalc(metaAMR90, amrT)

rgi7<-pvalCalc(metaRGI7, rgiT)
rgi10<-pvalCalc(metaRGI10, rgiT)
rgi30<-pvalCalc(metaRGI30, rgiT)
rgi60<-pvalCalc(metaRGI60, rgiT)
rgi90<-pvalCalc(metaRGI90, rgiT)

vsearch7<-pvalCalc(metaVSEARCH7, vsearchT)
vsearch10<-pvalCalc(metaVSEARCH10, vsearchT)
vsearch30<-pvalCalc(metaVSEARCH30, vsearchT)
vsearch60<-pvalCalc(metaVSEARCH60, vsearchT)
vsearch90<-pvalCalc(metaVSEARCH90, vsearchT)

#adjusting pvalues
bracken7adj<-apply(bracken7, 2, p.adjust, method = "BH")
bracken10adj<-apply(bracken10, 2, p.adjust, method = "BH")
bracken30adj<-apply(bracken30, 2, p.adjust, method = "BH")
bracken60adj<-apply(bracken60, 2, p.adjust, method = "BH")
bracken90adj<-apply(bracken90, 2, p.adjust, method = "BH")

amr7adj<-apply(amr7, 2, p.adjust, method = "BH")
amr10adj<-apply(amr10, 2, p.adjust, method = "BH")
amr30adj<-apply(amr30, 2, p.adjust, method = "BH")
amr60adj<-apply(amr60, 2, p.adjust, method = "BH")
amr90adj<-apply(amr90, 2, p.adjust, method = "BH")

rgi7adj<-apply(rgi7, 2, p.adjust, method = "BH")
rgi10adj<-apply(rgi10, 2, p.adjust, method = "BH")
rgi30adj<-apply(rgi30, 2, p.adjust, method = "BH")
rgi60adj<-apply(rgi60, 2, p.adjust, method = "BH")
rgi90adj<-apply(rgi90, 2, p.adjust, method = "BH")

vsearch7adj<-apply(vsearch7, 2, p.adjust, method = "BH")
vsearch10adj<-apply(vsearch10, 2, p.adjust, method = "BH")
vsearch30adj<-apply(vsearch30, 2, p.adjust, method = "BH")
vsearch60adj<-apply(vsearch60, 2, p.adjust, method = "BH")
vsearch90adj<-apply(vsearch90, 2, p.adjust, method = "BH")

write.csv(bracken7, "ABexposure/bracken7.csv")
write.csv(bracken10, "ABexposure/bracken10.csv")
write.csv(bracken30, "ABexposure/bracken30.csv")
write.csv(bracken60, "ABexposure/bracken60.csv")
write.csv(bracken90, "ABexposure/bracken90.csv")

write.csv(amr7, "ABexposure/amr7.csv")
write.csv(amr10, "ABexposure/amr10.csv")
write.csv(amr30, "ABexposure/amr30.csv")
write.csv(amr60, "ABexposure/amr60.csv")
write.csv(amr90, "ABexposure/amr90.csv")

write.csv(rgi7, "ABexposure/rgi7.csv")
write.csv(rgi10, "ABexposure/rgi10.csv")
write.csv(rgi30, "ABexposure/rgi30.csv")
write.csv(rgi60, "ABexposure/rgi60.csv")
write.csv(rgi90, "ABexposure/rgi90.csv")

write.csv(vsearch7, "ABexposure/vsearch7.csv")
write.csv(vsearch10, "ABexposure/vsearch10.csv")
write.csv(vsearch30, "ABexposure/vsearch30.csv")
write.csv(vsearch60, "ABexposure/vsearch60.csv")
write.csv(vsearch90, "ABexposure/vsearch90.csv")

write.csv(bracken7adj, "ABexposure/bracken7adj.csv")
write.csv(bracken10adj, "ABexposure/bracken10adj.csv")
write.csv(bracken30adj, "ABexposure/bracken30adj.csv")
write.csv(bracken60adj, "ABexposure/bracken60adj.csv")
write.csv(bracken90adj, "ABexposure/bracken90adj.csv")

write.csv(amr7adj, "ABexposure/amr7adj.csv")
write.csv(amr10adj, "ABexposure/amr10adj.csv")
write.csv(amr30adj, "ABexposure/amr30adj.csv")
write.csv(amr60adj, "ABexposure/amr60adj.csv")
write.csv(amr90adj, "ABexposure/amr90adj.csv")

write.csv(rgi7adj, "ABexposure/rgi7adj.csv")
write.csv(rgi10adj, "ABexposure/rgi10adj.csv")
write.csv(rgi30adj, "ABexposure/rgi30adj.csv")
write.csv(rgi60adj, "ABexposure/rgi60adj.csv")
write.csv(rgi90adj, "ABexposure/rgi90adj.csv")

write.csv(vsearch7adj, "ABexposure/vsearch7adj.csv")
write.csv(vsearch10adj, "ABexposure/vsearch10adj.csv")
write.csv(vsearch30adj, "ABexposure/vsearch30adj.csv")
write.csv(vsearch60adj, "ABexposure/vsearch60adj.csv")
write.csv(vsearch90adj, "ABexposure/vsearch90adj.csv")
