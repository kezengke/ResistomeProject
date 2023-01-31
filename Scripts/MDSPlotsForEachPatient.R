rm(list = ls())
library(stringr)
library("RColorBrewer")
library("vegan")
#metadata
dukeSamples<-read.csv("Duke_samples_meta.csv", header = T, sep = ",")
rownames(dukeSamples)<-dukeSamples$Seq_sample
rownames(dukeSamples)<-sapply(str_split(rownames(dukeSamples), "_", n = 2), `[`, 2)
dukeSamples$ID<-paste0("D", gsub("pre","" , sapply(str_split(dukeSamples$sample, "D", n = 3), `[`, 2), ignore.case = T)) #adding sample ID
dukeSamples$ID<-gsub("npe", "", dukeSamples$ID, ignore.case = T)

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

IDtypes<-unique(metaBracken$ID)
pdf("Plots/MDSPlotsForEachPatient(bracken).pdf", width=12, height=18)
par(mfrow=c(3,2))
par(mar=c(5,6,4,1)+.1)
for (i in 1:length(IDtypes)) {
  mdsT<-brackenT[, colnames(brackenT)[metaBracken$ID == IDtypes[i]]]
  if(ncol(mdsT)<3)
    next
  MDS<-capscale(t(mdsT)~1,distance = "bray")
  percentVariance<-MDS$CA$eig/sum(eigenvals(MDS))*100
  statusPlot<-ordiplot(MDS,choices = c(1,2),type="none",cex.lab = 1,
                       xlab = paste("MDS1  ", format(percentVariance[1], digits = 4), "%", sep = ""),
                       ylab = paste("MDS2  ", format(percentVariance[2], digits = 4), "%", sep = ""),
                       main = paste("Bracken-Species", IDtypes[i]))
  points(statusPlot,"sites", pch = 19, cex = 2.5, col = adjustcolor("tan2", alpha.f = 0.5))
}
dev.off()

IDtypes<-unique(metaAMR$ID)
pdf("Plots/MDSPlotsForEachPatient(AMR).pdf", width=12, height=18)
par(mfrow=c(3,2))
par(mar=c(5,6,4,1)+.1)
for (i in 1:length(IDtypes)) {
  mdsT<-amrT[, colnames(amrT)[metaAMR$ID == IDtypes[i]]]
  if(ncol(mdsT)<3)
    next
  MDS<-capscale(t(mdsT)~1,distance = "bray")
  percentVariance<-MDS$CA$eig/sum(eigenvals(MDS))*100
  statusPlot<-ordiplot(MDS,choices = c(1,2),type="none",cex.lab = 1,
                       xlab = paste("MDS1  ", format(percentVariance[1], digits = 4), "%", sep = ""),
                       ylab = paste("MDS2  ", format(percentVariance[2], digits = 4), "%", sep = ""),
                       main = paste("AMR", IDtypes[i]))
  points(statusPlot,"sites", pch = 19, cex = 2.5, col = adjustcolor("coral3", alpha.f = 0.5))
}
dev.off()

IDtypes<-unique(metaRGI$ID)
pdf("Plots/MDSPlotsForEachPatient(RGI).pdf", width=12, height=18)
par(mfrow=c(3,2))
par(mar=c(5,6,4,1)+.1)
for (i in 1:length(IDtypes)) {
  mdsT<-rgiT[, colnames(rgiT)[metaRGI$ID == IDtypes[i]]]
  if(ncol(mdsT)<3)
    next
  MDS<-capscale(t(mdsT)~1,distance = "bray")
  percentVariance<-MDS$CA$eig/sum(eigenvals(MDS))*100
  statusPlot<-ordiplot(MDS,choices = c(1,2),type="none",cex.lab = 1,
                       xlab = paste("MDS1  ", format(percentVariance[1], digits = 4), "%", sep = ""),
                       ylab = paste("MDS2  ", format(percentVariance[2], digits = 4), "%", sep = ""),
                       main = paste("RGI", IDtypes[i]))
  points(statusPlot,"sites", pch = 19, cex = 2.5, col = adjustcolor("cornflowerblue", alpha.f = 0.5))
}
dev.off()

IDtypes<-unique(metaVsearch$ID)
pdf("Plots/MDSPlotsForEachPatient(vsearch).pdf", width=12, height=18)
par(mfrow=c(3,2))
par(mar=c(5,6,4,1)+.1)
for (i in 1:length(IDtypes)) {
  mdsT<-vsearchT[, colnames(vsearchT)[metaVsearch$ID == IDtypes[i]]]
  if(ncol(mdsT)<3)
    next
  MDS<-capscale(t(mdsT)~1,distance = "bray")
  percentVariance<-MDS$CA$eig/sum(eigenvals(MDS))*100
  statusPlot<-ordiplot(MDS,choices = c(1,2),type="none",cex.lab = 1,
                       xlab = paste("MDS1  ", format(percentVariance[1], digits = 4), "%", sep = ""),
                       ylab = paste("MDS2  ", format(percentVariance[2], digits = 4), "%", sep = ""),
                       main = paste("vsearch", IDtypes[i]))
  points(statusPlot,"sites", pch = 19, cex = 2.5, col = adjustcolor("olivedrab4", alpha.f = 0.5))
}
dev.off()