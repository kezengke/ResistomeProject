rm(list = ls())
library(stringr)
library("RColorBrewer")
library("vegan")
library("dplyr")
#metadata
dukeSamples<-read.csv("Duke_samples_meta.csv", header = T, sep = ",")
rownames(dukeSamples)<-dukeSamples$Seq_sample
rownames(dukeSamples)<-sapply(str_split(rownames(dukeSamples), "_", n = 2), `[`, 2)
dukeSamples$ID<-paste0("D", gsub("pre","" , sapply(str_split(dukeSamples$sample, "D", n = 3), `[`, 2), ignore.case = T)) #adding sample ID
dukeSamples$ID<-gsub("npe", "", dukeSamples$ID, ignore.case = T)
dukeSamples<-dukeSamples %>% mutate(bins = case_when(between(dukeSamples$Timepoint, -100, -2) ~ "PRE",
                                                     # between(dukeSamples$Timepoint, -3, 3) ~ "D0",
                                                     between(dukeSamples$Timepoint, -3, 10) ~ "D7",
                                                     between(dukeSamples$Timepoint, 11, 17) ~ "D14",
                                                     between(dukeSamples$Timepoint, 18, 24) ~ "D21",
                                                     # between(dukeSamples$Timepoint, 25, 31) ~ "D28",
                                                     between(dukeSamples$Timepoint, 25, 45) ~ "D35",
                                                     between(dukeSamples$Timepoint, 46, 75) ~ "D60",
                                                     # between(dukeSamples$Timepoint, 76, 125) ~ "D100",
                                                     # between(dukeSamples$Timepoint, 126, 235) ~ "D180",
                                                     between(dukeSamples$Timepoint, 76, 965) ~ "D100"))

#bracken
brackenT<-read.delim("CountsTables/duke_bracken.csv", sep = ",", header = T, row.names = 1)
brackenT<-brackenT[, -c(1,2)] #get rid of taxonomy ID and level
brackenT<-brackenT[, grepl("num", colnames(brackenT))] #get rid of fractions

colnames(brackenT)<-gsub(".bracken.out_num", "", colnames(brackenT)) #get rid of useless info
SampleWithDots<-sapply(str_count(colnames(brackenT), "\\."), `[`, 1) == 2
colnames(brackenT)[SampleWithDots]<-sub("\\.", "", colnames(brackenT)[SampleWithDots]) #get rid of random "."s in sample names
colnames(brackenT)<-sub("\\.", "-", colnames(brackenT)) #replace . with - to match with metadata
colnames(brackenT)<-sapply(str_split(colnames(brackenT), "_", n = 2), `[`, 2) #keeping sequencing info for matching

#Normalization
n<-colSums(brackenT)
sumx<-sum(brackenT)
brackenT<-log10((brackenT/n)*(sumx/ncol(brackenT))+1)

#AMR
amrT<-read.delim("CountsTables/AMR_counts.tsv", sep = "\t", header = T, row.names = 1)
amrT<-amrT[, grepl("^D", colnames(amrT))]
colnames(amrT)<-gsub(".amrfinder.txt", "", colnames(amrT)) #get rid of useless info
SampleWithDots<-sapply(str_count(colnames(amrT), "\\."), `[`, 1) == 2
colnames(amrT)[SampleWithDots]<-sub("\\.", "", colnames(amrT)[SampleWithDots]) #get rid of random "."s in sample names
colnames(amrT)<-sub("\\.", "-", colnames(amrT)) #replace . with - to match with metadata
colnames(amrT)<-sapply(str_split(colnames(amrT), "_", n = 2), `[`, 2)

#RGI
rgiT<-read.delim("CountsTables/RGI_counts.tsv", sep = "\t", header = T, row.names = 1)
rgiT<-rgiT[, grepl("^D", colnames(rgiT))]
colnames(rgiT)<-gsub(".rgi.txt", "", colnames(rgiT)) #get rid of useless info
SampleWithDots<-sapply(str_count(colnames(rgiT), "\\."), `[`, 1) == 2
colnames(rgiT)[SampleWithDots]<-sub("\\.", "", colnames(rgiT)[SampleWithDots]) #get rid of random "."s in sample names
colnames(rgiT)<-sub("\\.", "-", colnames(rgiT)) #replace . with - to match with metadata
colnames(rgiT)<-sapply(str_split(colnames(rgiT), "_", n = 2), `[`, 2)

#vsearch
vsearchT<-read.delim("CountsTables/vsearch_counts.tsv", sep = "\t", header = T, row.names = 2)
vsearchT<-vsearchT[, -1] #get rid of index column
vsearchT<-vsearchT[, grepl("^D", colnames(vsearchT))]
colnames(vsearchT)<-gsub(".txt", "", colnames(vsearchT)) #get rid of useless info

SampleWithDots<-sapply(str_count(colnames(vsearchT), "\\."), `[`, 1) == 2
colnames(vsearchT)[SampleWithDots]<-sub("\\.", "", colnames(vsearchT)[SampleWithDots]) #get rid of random "."s in sample names
colnames(vsearchT)<-sub("\\.", "-", colnames(vsearchT)) #replace . with - to match with metadata
colnames(vsearchT)<-sapply(str_split(colnames(vsearchT), "_", n = 2), `[`, 2)

#meta for each table
metaBracken<-dukeSamples[colnames(brackenT), ]
metaAMR<-dukeSamples[colnames(amrT), ]
metaRGI<-dukeSamples[colnames(rgiT), ]
metaVsearch<-dukeSamples[colnames(vsearchT), ]

TimePointtypes<-unique(metaBracken$bins)
for (i in 1:length(TimePointtypes)) {
  mdsT<-brackenT[, colnames(brackenT)[metaBracken$bins == TimePointtypes[i]]]
  if(ncol(mdsT)<3)
    next
  MDS<-capscale(t(mdsT)~1,distance = "bray")
  percentVariance<-MDS$CA$eig/sum(eigenvals(MDS))*100
  pdf(paste("MDSPlotsForEachTimepoint/Bracken/Bracken-Species", TimePointtypes[i], ".pdf"), width=6, height=6)
  statusPlot<-ordiplot(MDS,choices = c(1,2),type="none",cex.lab = 1,
                       xlab = paste("MDS1  ", format(percentVariance[1], digits = 4), "%", sep = ""),
                       ylab = paste("MDS2  ", format(percentVariance[2], digits = 4), "%", sep = ""),
                       main = paste0("Bracken-Species ", TimePointtypes[i], " (n=", ncol(mdsT), ")"),
                       xlim = c(-2.5, 2.5), ylim = c(-2, 2))
  points(statusPlot,"sites", pch = 19, cex = 2.5, col = adjustcolor("tan2", alpha.f = 0.5))
  text(statusPlot, "sites", labels = dukeSamples[colnames(mdsT), 5], cex = 0.6, pos = 4, col = "tan2")
  dev.off()
}

TimePointtypes<-unique(metaAMR$bins)
for (i in 1:length(TimePointtypes)) {
  mdsT<-amrT[, colnames(amrT)[metaAMR$bins == TimePointtypes[i]]]
  if(ncol(mdsT)<3)
    next
  MDS<-capscale(t(mdsT)~1,distance = "bray")
  percentVariance<-MDS$CA$eig/sum(eigenvals(MDS))*100
  pdf(paste("MDSPlotsForEachTimepoint/AMR/AMR", TimePointtypes[i], ".pdf"), width=6, height=6)
  statusPlot<-ordiplot(MDS,choices = c(1,2),type="none",cex.lab = 1,
                       xlab = paste("MDS1  ", format(percentVariance[1], digits = 4), "%", sep = ""),
                       ylab = paste("MDS2  ", format(percentVariance[2], digits = 4), "%", sep = ""),
                       main = paste0("AMR ", TimePointtypes[i], " (n=", ncol(mdsT), ")"),
                       xlim = c(-3, 2.5), ylim = c(-2, 2.5))
  points(statusPlot,"sites", pch = 19, cex = 2.5, col = adjustcolor("coral3", alpha.f = 0.5))
  text(statusPlot, "sites", labels = dukeSamples[colnames(mdsT), 5], cex = 0.6, pos = 4, col = "coral3")
  dev.off()
}

TimePointtypes<-unique(metaRGI$bins)
for (i in 1:length(TimePointtypes)) {
  mdsT<-rgiT[, colnames(rgiT)[metaRGI$bins == TimePointtypes[i]]]
  if(ncol(mdsT)<3)
    next
  MDS<-capscale(t(mdsT)~1,distance = "bray")
  percentVariance<-MDS$CA$eig/sum(eigenvals(MDS))*100
  pdf(paste("MDSPlotsForEachTimepoint/RGI/RGI", TimePointtypes[i], ".pdf"), width=6, height=6)
  statusPlot<-ordiplot(MDS,choices = c(1,2),type="none",cex.lab = 1,
                       xlab = paste("MDS1  ", format(percentVariance[1], digits = 4), "%", sep = ""),
                       ylab = paste("MDS2  ", format(percentVariance[2], digits = 4), "%", sep = ""),
                       main = paste0("RGI ", TimePointtypes[i], " (n=", ncol(mdsT), ")"),
                       xlim = c(-2, 2), ylim = c(-2, 3))
  points(statusPlot,"sites", pch = 19, cex = 2.5, col = adjustcolor("cornflowerblue", alpha.f = 0.5))
  text(statusPlot, "sites", labels = dukeSamples[colnames(mdsT), 5], cex = 0.6, pos = 4, col = "cornflowerblue")
  dev.off()
}

TimePointtypes<-unique(metaVsearch$bins)
for (i in 1:length(TimePointtypes)) {
  mdsT<-vsearchT[, colnames(vsearchT)[metaVsearch$bins == TimePointtypes[i]]]
  if(ncol(mdsT)<3)
    next
  MDS<-capscale(t(mdsT)~1,distance = "bray")
  percentVariance<-MDS$CA$eig/sum(eigenvals(MDS))*100
  pdf(paste("MDSPlotsForEachTimepoint/vsearch/vsearch", TimePointtypes[i], ".pdf"), width=6, height=6)
  statusPlot<-ordiplot(MDS,choices = c(1,2),type="none",cex.lab = 1,
                       xlab = paste("MDS1  ", format(percentVariance[1], digits = 4), "%", sep = ""),
                       ylab = paste("MDS2  ", format(percentVariance[2], digits = 4), "%", sep = ""),
                       main = paste0("vsearch ", TimePointtypes[i], " (n=", ncol(mdsT), ")"),
                       xlim = c(-1.5, 2.5), ylim = c(-2, 2))
  points(statusPlot,"sites", pch = 19, cex = 2.5, col = adjustcolor("olivedrab4", alpha.f = 0.5))
  text(statusPlot, "sites", labels = dukeSamples[colnames(mdsT), 5], cex = 0.6, pos = 4, col = "olivedrab4")
  dev.off()
}
