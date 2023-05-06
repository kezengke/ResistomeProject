rm(list = ls())
library(stringr)
library("RColorBrewer")
library("vegan")
library("dplyr")
#metadata
dukeSamples<-read.csv("Duke_samples_meta.csv", header = T, sep = ",")
rownames(dukeSamples)<-dukeSamples$Seq_sample
rownames(dukeSamples)<-sapply(str_split(rownames(dukeSamples), "_", n = 2), `[`, 2)

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
for (i in 1:ncol(brackenT)) {
        brackenT[,i]<-brackenT[,i]/n[i]
}
brackenT<-log10(brackenT*(sumx/ncol(brackenT))+1)

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

#Bracken
pdf("Plots/TimeVsMDS1234(Bracken).pdf", width=12, height=12)
par(mfrow=c(2,2))
par(mar=c(5,6,4,1)+.1)
MDS<-summary(capscale(t(brackenT)~1,distance = "bray"))[["sites"]][,c(1,2,3,4)]
MDS<-log10(MDS+10)
plot(dukeSamples[rownames(MDS), 2], MDS[ ,1], 
     col = "tan2", pch = 19, 
     xlab = "Days", ylab = "log10(MDS1+10)",
     main = "Bracken MDS1")

plot(dukeSamples[rownames(MDS), 2], MDS[ ,2], 
     col = "tan2", pch = 19, 
     xlab = "Days", ylab = "log10(MDS2+10)",
     main = "Bracken MDS2")

plot(dukeSamples[rownames(MDS), 2], MDS[ ,3], 
     col = "tan2", pch = 19, 
     xlab = "Days", ylab = "log10(MDS3+10)",
     main = "Bracken MDS3")

plot(dukeSamples[rownames(MDS), 2], MDS[ ,4], 
     col = "tan2", pch = 19, 
     xlab = "Days", ylab = "log10(MDS4+10)",
     main = "Bracken MDS4")
dev.off()

#AMR
pdf("Plots/TimeVsMDS1234(AMR).pdf", width=12, height=12)
par(mfrow=c(2,2))
par(mar=c(5,6,4,1)+.1)
MDS<-summary(capscale(t(amrT)~1,distance = "bray"))[["sites"]][,c(1,2,3,4)]
MDS<-log10(MDS+10)
plot(dukeSamples[rownames(MDS), 2], MDS[ ,1], 
     col = "coral3", pch = 19, 
     xlab = "Days", ylab = "log10(MDS1+10)",
     main = "AMR MDS1")

plot(dukeSamples[rownames(MDS), 2], MDS[ ,2], 
     col = "coral3", pch = 19, 
     xlab = "Days", ylab = "log10(MDS2+10)",
     main = "AMR MDS2")

plot(dukeSamples[rownames(MDS), 2], MDS[ ,3], 
     col = "coral3", pch = 19, 
     xlab = "Days", ylab = "log10(MDS3+10)",
     main = "AMR MDS3")

plot(dukeSamples[rownames(MDS), 2], MDS[ ,4], 
     col = "coral3", pch = 19, 
     xlab = "Days", ylab = "log10(MDS4+10)",
     main = "AMR MDS4")
dev.off()

#AMR
pdf("Plots/TimeVsMDS1234(RGI).pdf", width=12, height=12)
par(mfrow=c(2,2))
par(mar=c(5,6,4,1)+.1)
MDS<-summary(capscale(t(rgiT)~1,distance = "bray"))[["sites"]][,c(1,2,3,4)]
MDS<-log10(MDS+10)
plot(dukeSamples[rownames(MDS), 2], MDS[ ,1], 
     col = "cornflowerblue", pch = 19, 
     xlab = "Days", ylab = "log10(MDS1+10)",
     main = "RGI MDS1")

plot(dukeSamples[rownames(MDS), 2], MDS[ ,2], 
     col = "cornflowerblue", pch = 19, 
     xlab = "Days", ylab = "log10(MDS2+10)",
     main = "RGI MDS2")

plot(dukeSamples[rownames(MDS), 2], MDS[ ,3], 
     col = "cornflowerblue", pch = 19, 
     xlab = "Days", ylab = "log10(MDS3+10)",
     main = "RGI MDS3")

plot(dukeSamples[rownames(MDS), 2], MDS[ ,4], 
     col = "cornflowerblue", pch = 19, 
     xlab = "Days", ylab = "log10(MDS4+10)",
     main = "RGI MDS4")
dev.off()

#vsearch
pdf("Plots/TimeVsMDS1234(vsearch).pdf", width=12, height=12)
par(mfrow=c(2,2))
par(mar=c(5,6,4,1)+.1)
MDS<-summary(capscale(t(vsearchT)~1,distance = "bray"))[["sites"]][,c(1,2,3,4)]
MDS<-log10(MDS+10)
plot(dukeSamples[rownames(MDS), 2], MDS[ ,1], 
     col = "olivedrab4", pch = 19, 
     xlab = "Days", ylab = "log10(MDS1+10)",
     main = "vsearch MDS1")

plot(dukeSamples[rownames(MDS), 2], MDS[ ,2], 
     col = "olivedrab4", pch = 19, 
     xlab = "Days", ylab = "log10(MDS2+10)",
     main = "vsearch MDS2")

plot(dukeSamples[rownames(MDS), 2], MDS[ ,3], 
     col = "olivedrab4", pch = 19, 
     xlab = "Days", ylab = "log10(MDS3+10)",
     main = "vsearch MDS3")

plot(dukeSamples[rownames(MDS), 2], MDS[ ,4], 
     col = "olivedrab4", pch = 19, 
     xlab = "Days", ylab = "log10(MDS4+10)",
     main = "vsearch MDS4")
dev.off()

