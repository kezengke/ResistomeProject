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

pdf("Plots/TotalCountsByTimePointsBins.pdf", width=12, height=12)
par(mfrow=c(2, 2))
par(mar=c(5, 6, 4, 1)+.1)

binCounts<-aggregate(colSums(brackenT), list(dukeSamples[colnames(brackenT), 5]), FUN=sum)
binCounts<-binCounts[c(11, 1, 10, 3, 5, 6, 7, 9, 2, 4, 8),]
barplot(binCounts$x, names.arg=c("PRE", "D0", "D7", "D14","D21" , "D28", "D35", "D60", "D100", "D180", "D365"), 
        main = "Taxa Counts by Time (Bracken-Species)", cex.axis = 1.5, cex.lab = 1.5, cex.names = 0.8, cex.main = 1.8,
        xlab = "Time Points (Days)", ylab = "Counts(Normalized)", col = "tan2")

binCounts<-aggregate(colSums(amrT), list(dukeSamples[colnames(amrT), 5]), FUN=sum)
binCounts<-binCounts[c(11, 1, 10, 3, 5, 6, 7, 9, 2, 4, 8),]
barplot(binCounts$x, names.arg=c("PRE", "D0", "D7", "D14","D21" , "D28", "D35", "D60", "D100", "D180", "D365"), 
        main = "Gene Counts by Time (AMR)", cex.axis = 1.5, cex.lab = 1.5, cex.names = 0.8, cex.main = 1.8,
        xlab = "Time Points (Days)", ylab = "Counts", col = "coral3")

binCounts<-aggregate(colSums(rgiT), list(dukeSamples[colnames(rgiT), 5]), FUN=sum)
binCounts<-binCounts[c(11, 1, 10, 3, 5, 6, 7, 9, 2, 4, 8),]
barplot(binCounts$x, names.arg=c("PRE", "D0", "D7", "D14","D21" , "D28", "D35", "D60", "D100", "D180", "D365"), 
        main = "Gene Counts by Time (RGI)", cex.axis = 1.5, cex.lab = 1.5, cex.names = 0.8, cex.main = 1.8,
        xlab = "Time Points (Days)", ylab = "Counts", col = "cornflowerblue")

binCounts<-aggregate(colSums(vsearchT), list(dukeSamples[colnames(vsearchT), 5]), FUN=sum)
binCounts<-binCounts[c(11, 1, 10, 3, 5, 6, 7, 9, 2, 4, 8),]
barplot(binCounts$x, names.arg=c("PRE", "D0", "D7", "D14","D21" , "D28", "D35", "D60", "D100", "D180", "D365"), 
        main = "Gene Counts by Time (vsearch)", cex.axis = 1.5, cex.lab = 1.5, cex.names = 0.8, cex.main = 1.8,
        xlab = "Time Points (Days)", ylab = "Counts", col = "olivedrab4")

dev.off()
