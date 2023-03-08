rm(list = ls())
library(stringr)
library("RColorBrewer")
library("vegan")
library("dplyr")
library("nlme")
#metadata
dukeSamples<-read.csv("Duke_samples_meta.csv", header = T, sep = ",")
rownames(dukeSamples)<-dukeSamples$Seq_sample
rownames(dukeSamples)<-sapply(str_split(rownames(dukeSamples), "_", n = 2), `[`, 2)
dukeSamples$ID<-paste0("D", gsub("pre","" , sapply(str_split(dukeSamples$sample, "D", n = 3), `[`, 2), ignore.case = T)) #adding sample ID
dukeSamples$ID<-gsub("npe", "", dukeSamples$ID, ignore.case = T)
dukeSamples<-dukeSamples %>% mutate(bins = case_when(between(dukeSamples$Timepoint, -100, -2) ~ "PRE",
                                                     between(dukeSamples$Timepoint, -3, 10) ~ "D7",
                                                     between(dukeSamples$Timepoint, 11, 17) ~ "D14",
                                                     between(dukeSamples$Timepoint, 18, 24) ~ "D21",
                                                     between(dukeSamples$Timepoint, 25, 45) ~ "D35",
                                                     between(dukeSamples$Timepoint, 46, 75) ~ "D60",
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

#AMR
amrT<-read.delim("CountsTables/AMR_counts.tsv", sep = "\t", header = T, row.names = 2)
amrT<-amrT[,-1]
amrT<-amrT[, grepl("^D", colnames(amrT))]
colnames(amrT)<-gsub(".amrfinder.txt", "", colnames(amrT)) #get rid of useless info
SampleWithDots<-sapply(str_count(colnames(amrT), "\\."), `[`, 1) == 2
colnames(amrT)[SampleWithDots]<-sub("\\.", "", colnames(amrT)[SampleWithDots]) #get rid of random "."s in sample names
colnames(amrT)<-sub("\\.", "-", colnames(amrT)) #replace . with - to match with metadata
colnames(amrT)<-sapply(str_split(colnames(amrT), "_", n = 2), `[`, 2)

#RGI
rgiT<-read.delim("CountsTables/RGI_counts.tsv", sep = "\t", header = T, row.names = 2)
rgiT<-rgiT[, -1]
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


#Normalization
n<-colSums(brackenT)
sumx<-sum(brackenT)
brackenT<-log10((brackenT/n)*(sumx/ncol(brackenT))+1)

#lm 2nd order
pdf("Plots/SecondOrderLm(Bracken).pdf", width=12, height=18)
par(mfrow=c(3, 2))
par(mar=c(5, 6, 4, 1)+.1)
brackenT<-brackenT[apply(brackenT == 0, 1, sum) <= (ncol(brackenT)*0.8), ]
for (i in 1:nrow(brackenT)) {
  Model<-lm(unlist(brackenT[i,]) ~ poly(dukeSamples$Timepoint, 2) , x=T)
  plot(dukeSamples$Timepoint, unlist(brackenT[i,]), 
       pch = 19, col = "tan2", 
       main = paste("Bracken(Species)\n", rownames(brackenT)[i]),
       xlab = "TimePoints", ylab = "Counts(Log10)")
  lines(sort(dukeSamples$Timepoint), fitted(Model)[order(dukeSamples$Timepoint)], 
        col = "grey", type = "l")
}
dev.off()

pdf("Plots/SecondOrderLm(AMR).pdf", width=12, height=18)
par(mfrow=c(3, 2))
par(mar=c(5, 6, 4, 1)+.1)
amrT<-amrT[apply(amrT == 0, 1, sum) <= (ncol(amrT)*0.8), ]
for (i in 1:nrow(amrT)) {
  Model<-lm(unlist(amrT[i,]) ~ poly(dukeSamples$Timepoint, 2) , x=T)
  plot(dukeSamples$Timepoint, unlist(amrT[i,]), 
       pch = 19, col = "coral3", 
       main = paste("AMR\n", rownames(amrT)[i]),
       xlab = "TimePoints", ylab = "Counts")
  lines(sort(dukeSamples$Timepoint), fitted(Model)[order(dukeSamples$Timepoint)], 
        col = "grey", type = "l")
}
dev.off()

pdf("Plots/SecondOrderLm(RGI).pdf", width=12, height=18)
par(mfrow=c(3, 2))
par(mar=c(5, 6, 4, 1)+.1)
rgiT<-rgiT[apply(rgiT == 0, 1, sum) <= (ncol(rgiT)*0.8), ]
for (i in 1:nrow(rgiT)) {
  Model<-lm(unlist(rgiT[i,]) ~ poly(dukeSamples$Timepoint, 2) , x=T)
  plot(dukeSamples$Timepoint, unlist(rgiT[i,]), 
       pch = 19, col = "cornflowerblue", 
       main = paste("RGI\n", rownames(rgiT)[i]),
       xlab = "TimePoints", ylab = "Counts")
  lines(sort(dukeSamples$Timepoint), fitted(Model)[order(dukeSamples$Timepoint)], 
        col = "grey", type = "l")
}
dev.off()

pdf("Plots/SecondOrderLm(vsearch).pdf", width=12, height=18)
par(mfrow=c(3, 2))
par(mar=c(5, 6, 4, 1)+.1)
vsearchT<-vsearchT[apply(vsearchT == 0, 1, sum) <= (ncol(vsearchT)*0.8), ]
for (i in 1:nrow(vsearchT)) {
  Model<-lm(unlist(vsearchT[i,]) ~ poly(dukeSamples$Timepoint, 2) , x=T)
  plot(dukeSamples$Timepoint, unlist(vsearchT[i,]), 
       pch = 19, col = "olivedrab4", 
       main = paste("vsearch\n", rownames(vsearchT)[i]),
       xlab = "TimePoints", ylab = "Counts")
  lines(sort(dukeSamples$Timepoint), fitted(Model)[order(dukeSamples$Timepoint)], 
        col = "grey", type = "l")
}
dev.off()