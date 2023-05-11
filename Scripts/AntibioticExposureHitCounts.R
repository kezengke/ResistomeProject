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

#bracken
brackenT<-read.csv("CountsTables/brackenNormalized.csv", header = T, row.names = 1, check.names = F)
amrT<-read.csv("CountsTables/amrNormalized.csv", header = T, row.names = 1, check.names = F)
rgiT<-read.csv("CountsTables/rgiNormalized.csv", header = T, row.names = 1, check.names = F)
vsearchT<-read.csv("CountsTables/vsearchNormalized.csv", header = T, row.names = 1, check.names = F)
#filter
brackenT<-brackenT[apply(brackenT == 0, 1, sum) <= (ncol(brackenT)*0.8), ]
amrT<-amrT[apply(amrT == 0, 1, sum) <= (ncol(amrT)*0.8), ]
rgiT<-rgiT[apply(rgiT == 0, 1, sum) <= (ncol(rgiT)*0.8), ]
vsearchT<-vsearchT[apply(vsearchT == 0, 1, sum) <= (ncol(vsearchT)*0.8), ]

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

hitsCalc<-function(meta, countsT) {
  hits<-vector()
  for (i in 1:ncol(meta)) {
    testingMeta<-meta[, i, drop = F]
    pval<-vector()
    for (j in 1:nrow(countsT)) {
      if (colSums(testingMeta == "Yes") < 2) {
        pval[j]<-NA
      }
      else
        pval[j]<-t.test(unlist(countsT[j, ])~factor(testingMeta[,1]))$p.value
    }
    hits[i]<-sum(pval<0.05)
  }
  return(hits)
}

bracken7<-hitsCalc(metaBRACKEN7, brackenT)
bracken10<-hitsCalc(metaBRACKEN10, brackenT)
bracken30<-hitsCalc(metaBRACKEN30, brackenT)
bracken60<-hitsCalc(metaBRACKEN60, brackenT)
bracken90<-hitsCalc(metaBRACKEN90, brackenT)

amr7<-hitsCalc(metaAMR7, amrT)
amr10<-hitsCalc(metaAMR10, amrT)
amr30<-hitsCalc(metaAMR30, amrT)
amr60<-hitsCalc(metaAMR60, amrT)
amr90<-hitsCalc(metaAMR90, amrT)

rgi7<-hitsCalc(metaRGI7, rgiT)
rgi10<-hitsCalc(metaRGI10, rgiT)
rgi30<-hitsCalc(metaRGI30, rgiT)
rgi60<-hitsCalc(metaRGI60, rgiT)
rgi90<-hitsCalc(metaRGI90, rgiT)

vsearch7<-hitsCalc(metaVSEARCH7, vsearchT)
vsearch10<-hitsCalc(metaVSEARCH10, vsearchT)
vsearch30<-hitsCalc(metaVSEARCH30, vsearchT)
vsearch60<-hitsCalc(metaVSEARCH60, vsearchT)
vsearch90<-hitsCalc(metaVSEARCH90, vsearchT)

HitsT<-rbind(bracken7, bracken10, bracken30, bracken60, bracken90, 
             amr7, amr10, amr30, amr60, amr90, 
             rgi7, rgi10, rgi30, rgi60, rgi90,
             vsearch7, vsearch10, vsearch30, vsearch60, vsearch90)

colnames(HitsT)<-colnames(metaBRACKEN7)

write.csv(HitsT, "AntibioticExposureHitCounts.csv")
