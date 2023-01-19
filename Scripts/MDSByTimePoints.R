rm(list = ls())
library(stringr)
library("RColorBrewer")
library("vegan")
library("dplyr")

#metadata
dukeSamples<-read.csv("Duke_samples_meta.csv", header = T, sep = ",")
rownames(dukeSamples)<-dukeSamples$Seq_sample
rownames(dukeSamples)<-sapply(str_split(rownames(dukeSamples), "_", n = 2), `[`, 2)

dukeSamples<-dukeSamples %>% mutate(bins = case_when(dukeSamples$Timepoint < 1 ~ "PRE",
                                        between(dukeSamples$Timepoint, 1, 7) ~ "Week1",
                                        between(dukeSamples$Timepoint, 8, 14) ~ "Week2",
                                        between(dukeSamples$Timepoint, 15, 21) ~ "Week3",
                                        between(dukeSamples$Timepoint, 22, 28) ~ "Week4",
                                        dukeSamples$Timepoint > 28 ~ "AfterDay28"))

#bracken
brackenT<-read.delim("CountsTables/duke_bracken.csv", sep = ",", header = T, row.names = 1)
brackenT<-brackenT[, -c(1,2)] #get rid of taxonomy ID and level
brackenT<-brackenT[, grepl("num", colnames(brackenT))] #get rid of fractions

colnames(brackenT)<-gsub(".bracken.out_num", "", colnames(brackenT)) #get rid of useless info
SampleWithDots<-sapply(str_count(colnames(brackenT), "\\."), `[`, 1) == 2
colnames(brackenT)[SampleWithDots]<-sub("\\.", "", colnames(brackenT)[SampleWithDots]) #get rid of random "."s in sample names
colnames(brackenT)<-sub("\\.", "-", colnames(brackenT)) #replace . with - to match with metadata
colnames(brackenT)<-sapply(str_split(colnames(brackenT), "_", n = 2), `[`, 2) #keeping sequencing info for matching

dukeSamples[setdiff(rownames(dukeSamples), colnames(brackenT)), 4] #for missing samples in bracken counts table 

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

pdf("Plots/MDSForAll(ColoredByTimePoints).pdf", width=12, height=12)
par(mfrow=c(2,2))
par(mar=c(5,6,4,1)+.1)

#Bracken
circleCol<-brewer.pal(length(unique(metaBracken$bins)), "Set3")
cols<-circleCol[factor(metaBracken$bins)]

MDS<-capscale(t(brackenT)~1,distance = "bray")
percentVariance<-MDS$CA$eig/sum(eigenvals(MDS))*100
cols<-circleCol
statusPlot<-ordiplot(MDS,choices = c(1,2),type="none",cex.lab = 1,
                     xlab = paste("MDS1  ", format(percentVariance[1], digits = 4), "%", sep = ""),
                     ylab = paste("MDS2  ", format(percentVariance[2], digits = 4), "%", sep = ""),
                     main = paste("Bracken-Species n =", ncol(brackenT)))
points(statusPlot,"sites", pch = 19, cex = 2.5, col = adjustcolor(cols, alpha.f = 0.5))
ordiellipse(statusPlot, metaBracken$bins, kind="se", conf=0.95, lwd=4, draw = "lines", col=circleCol) 
legend("topright", sort(unique(metaBracken$bins)), col = circleCol[1:6], cex = 1, pch = 16, bty = "n")

#AMR
circleCol<-brewer.pal(length(unique(metaAMR$bins)), "Set3")
cols<-circleCol[factor(metaAMR$bins)]

MDS<-capscale(t(amrT)~1,distance = "bray")
percentVariance<-MDS$CA$eig/sum(eigenvals(MDS))*100
cols<-circleCol
statusPlot<-ordiplot(MDS,choices=c(1,2),type="none",cex.lab=1,
                     xlab=paste("MDS1  ", format(percentVariance[1], digits = 4), "%", sep = ""),
                     ylab=paste("MDS2  ", format(percentVariance[2], digits = 4), "%", sep = ""),
                     main=paste("AMR n =", ncol(amrT)))
points(statusPlot,"sites", pch=19, cex=2.5, col=adjustcolor(cols, alpha.f = 0.5))
ordiellipse(statusPlot, metaAMR$bins, kind="se", conf=0.95, lwd=4, draw = "lines", col=circleCol) 
legend("topleft", sort(unique(metaAMR$bins)), col = circleCol[1:6], cex = 1, pch = 16, bty = "n")

#RGI
circleCol<-brewer.pal(length(unique(metaRGI$bins)), "Set3")
cols<-circleCol[factor(metaRGI$bins)]

MDS<-capscale(t(rgiT)~1,distance = "bray")
percentVariance<-MDS$CA$eig/sum(eigenvals(MDS))*100
cols<-circleCol
statusPlot<-ordiplot(MDS,choices=c(1,2),type="none",cex.lab=1,
                     xlab=paste("MDS1  ", format(percentVariance[1], digits = 4), "%", sep = ""),
                     ylab=paste("MDS2  ", format(percentVariance[2], digits = 4), "%", sep = ""),
                     main=paste("RGI n =", ncol(rgiT)))
points(statusPlot,"sites", pch=19, cex=2.5, col=adjustcolor(cols, alpha.f = 0.5))
ordiellipse(statusPlot, metaRGI$bins, kind="se", conf=0.95, lwd=4, draw = "lines", col=circleCol) 
legend("topright", sort(unique(metaRGI$bins)), col = circleCol[1:6], cex = 1, pch = 16, bty = "n")

#vsearch
circleCol<-brewer.pal(length(unique(metaVsearch$bins)), "Set3") 
cols<-circleCol[factor(metaVsearch$bins)]

MDS<-capscale(t(vsearchT)~1,distance = "bray")
percentVariance<-MDS$CA$eig/sum(eigenvals(MDS))*100
cols<-circleCol
statusPlot<-ordiplot(MDS,choices=c(1,2),type="none",cex.lab=1,
                     xlab=paste("MDS1  ", format(percentVariance[1], digits = 4), "%", sep = ""),
                     ylab=paste("MDS2  ", format(percentVariance[2], digits = 4), "%", sep = ""),
                     main=paste("vsearch n =", ncol(vsearchT)))
points(statusPlot,"sites", pch=19, cex=2.5, col=adjustcolor(cols, alpha.f = 0.5))
ordiellipse(statusPlot, metaVsearch$bins, kind="se", conf=0.95, lwd=4, draw = "lines", col=circleCol) 
legend("topright", sort(unique(metaVsearch$bins)), col = circleCol[1:6], cex = 1, pch = 16, bty = "n")

dev.off()