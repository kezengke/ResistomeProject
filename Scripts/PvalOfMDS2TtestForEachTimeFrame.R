#-log10 pval over time for each timeframe MDS2
rm(list = ls())
library("vegan")

#metadata
metaData<-read.csv("metaWithBins.csv", header = T, row.names = 1)

#gene counts tables
speciesT<-read.csv("CountsTables/brackenFiltered.csv", header = T, row.names = 1, check.names = F)
genusT<-read.csv("CountsTables/genusFiltered.csv", header = T, row.names = 1, check.names = F)
amrT<-read.csv("CountsTables/amrFiltered.csv", header = T, row.names = 1, check.names = F)
rgiT<-read.csv("CountsTables/rgiFiltered.csv", header = T, row.names = 1, check.names = F)
vsearchT<-read.csv("CountsTables/vsearchFiltered.csv", header = T, row.names = 1, check.names = F)

timeFrames<-unique(metaData$bins)
desired_order<-c("PRE", "D0", "D7", "D14", "D21", "D28", "D35", "D60", "D100", "D180")
ordered_string<-factor(timeFrames, levels = desired_order, ordered = TRUE)
reordered_string<-timeFrames[order(ordered_string)]
timeFrames<-reordered_string

pvalspecies<-vector()
pvalgenus<-vector()
pvalamr<-vector()
pvalrgi<-vector()
pvalvsearch<-vector()
for (i in 1:length(timeFrames)) {
  samples<-rownames(metaData)[metaData$bins == timeFrames[i]]
  species<-speciesT[, samples, drop = F]
  genus<-genusT[, samples, drop = F]
  amr<-amrT[, intersect(colnames(amrT), samples), drop = F]
  rgi<-rgiT[, samples, drop = F]
  vsearch<-vsearchT[, samples, drop = F]
  
  metaSPECIES<-metaData[colnames(species), , drop = F]
  metaGENUS<-metaData[colnames(genus), , drop = F]
  metaAMR<-metaData[colnames(amr), , drop = F]
  metaRGI<-metaData[colnames(rgi), , drop = F]
  metaVSEARCH<-metaData[colnames(vsearch), , drop = F]
  
  MDSspecies<-summary(capscale(t(species)~1,distance = "bray"))[["sites"]][,2]
  MDSgenus<-summary(capscale(t(genus)~1,distance = "bray"))[["sites"]][,2]
  MDSamr<-summary(capscale(t(amr)~1,distance = "bray"))[["sites"]][,2]
  MDSrgi<-summary(capscale(t(rgi)~1,distance = "bray"))[["sites"]][,2]
  MDSvsearch<-summary(capscale(t(vsearch)~1,distance = "bray"))[["sites"]][,2]
  
  pvalspecies[i]<-log10(t.test(MDSspecies~metaSPECIES$ptInOut)$p.value) * (-1)
  pvalgenus[i]<-log10(t.test(MDSgenus~metaGENUS$ptInOut)$p.value) * (-1)
  pvalamr[i]<-log10(t.test(MDSamr~metaAMR$ptInOut)$p.value) * (-1)
  pvalrgi[i]<-log10(t.test(MDSrgi~metaRGI$ptInOut)$p.value) * (-1)
  pvalvsearch[i]<-log10(t.test(MDSvsearch~metaVSEARCH$ptInOut)$p.value) * (-1)
  
} 

pdf("Plots/TtestPvalOfMDS2ForEachTimeFramePlots.pdf", width=15, height=10)
par(mfrow=c(2, 3))
par(mar=c(5, 6, 4, 1)+.1)

plot(pvalspecies, xlab = "TimeFrame", ylab = "-log10(pvalue) of MDS2 t-test",
     xaxt = "n", col = "tan2", pch = 19, main = "Species(In vs. Out Patient)")
axis(1, at = 1:length(desired_order), labels = desired_order)

plot(pvalgenus, xlab = "TimeFrame", ylab = "-log10(pvalue) of MDS2 t-test",
     xaxt = "n", col = "orchid3", pch = 19, main = "Genus(In vs. Out Patient)")
axis(1, at = 1:length(desired_order), labels = desired_order)

plot(pvalamr, xlab = "TimeFrame", ylab = "-log10(pvalue) of MDS2 t-test",
     xaxt = "n", col = "coral3", pch = 19, main = "AMR(In vs. Out Patient)")
axis(1, at = 1:length(desired_order), labels = desired_order)

plot(pvalrgi, xlab = "TimeFrame", ylab = "-log10(pvalue) of MDS2 t-test",
     xaxt = "n", col = "cornflowerblue", pch = 19, main = "RGI(In vs. Out Patient)")
axis(1, at = 1:length(desired_order), labels = desired_order)

plot(pvalvsearch, xlab = "TimeFrame", ylab = "-log10(pvalue) of MDS2 t-test",
     xaxt = "n", col = "olivedrab4", pch = 19, main = "Vsearch(In vs. Out Patient)")
axis(1, at = 1:length(desired_order), labels = desired_order)

dev.off()
