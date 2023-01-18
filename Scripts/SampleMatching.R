rm(list = ls())
library(stringr)
dukeSamples<-read.csv("Duke_samples_meta.csv", header = T, sep = ",")
dukeSamples$sample<-gsub("pre", "PRE", dukeSamples$sample, ignore.case = T) #switch all pre to cap
dukeSamples$Seq_sample<-gsub("pre", "PRE", dukeSamples$Seq_sample, ignore.case = T) #switch all pre to cap

#AMR
amrT<-read.delim("CountsTables/AMR_counts.tsv", sep = "\t", header = T, row.names = 1)
amrT<-amrT[, grepl("^D", colnames(amrT))]
colnames(amrT)<-gsub("pre", "PRE", colnames(amrT), ignore.case = T) #switch all pre to cap
colnames(amrT)<-gsub(".amrfinder.txt", "", colnames(amrT)) #get rid of useless info

SampleWithDots<-sapply(str_count(colnames(amrT), "\\."), `[`, 1) == 2
colnames(amrT)[SampleWithDots]<-sub("\\.", "", colnames(amrT)[SampleWithDots]) #get rid of random "."s in sample names

colnames(amrT)<-sub("\\.", "-", colnames(amrT)) #replace . with - to match with metadata

length(intersect(dukeSamples$Seq_sample, colnames(amrT)))

setdiff(dukeSamples$Seq_sample, colnames(amrT))
setdiff(colnames(amrT), dukeSamples$Seq_sample)

#RGI
rgiT<-read.delim("CountsTables/RGI_counts.tsv", sep = "\t", header = T, row.names = 1)
rgiT<-rgiT[, grepl("^D", colnames(rgiT))]
colnames(rgiT)<-gsub("pre", "PRE", colnames(rgiT), ignore.case = T) #switch all pre to cap
colnames(rgiT)<-gsub(".rgi.txt", "", colnames(rgiT)) #get rid of useless info

SampleWithDots<-sapply(str_count(colnames(rgiT), "\\."), `[`, 1) == 2
colnames(rgiT)[SampleWithDots]<-sub("\\.", "", colnames(rgiT)[SampleWithDots]) #get rid of random "."s in sample names

colnames(rgiT)<-sub("\\.", "-", colnames(rgiT)) #replace . with - to match with metadata

length(intersect(dukeSamples$Seq_sample, colnames(rgiT)))

setdiff(dukeSamples$Seq_sample, colnames(rgiT))
setdiff(colnames(rgiT), dukeSamples$Seq_sample)

#vsearch
vsearchT<-read.delim("CountsTables/vsearch_counts.tsv", sep = "\t", header = T, row.names = 1)
vsearchT<-vsearchT[, grepl("^D", colnames(vsearchT))]
colnames(vsearchT)<-gsub("pre", "PRE", colnames(vsearchT), ignore.case = T) #switch all pre to cap
colnames(vsearchT)<-gsub(".txt", "", colnames(vsearchT)) #get rid of useless info

SampleWithDots<-sapply(str_count(colnames(vsearchT), "\\."), `[`, 1) == 2
colnames(vsearchT)[SampleWithDots]<-sub("\\.", "", colnames(vsearchT)[SampleWithDots]) #get rid of random "."s in sample names

colnames(vsearchT)<-sub("\\.", "-", colnames(vsearchT)) #replace . with - to match with metadata

length(intersect(dukeSamples$Seq_sample, colnames(vsearchT)))

setdiff(dukeSamples$Seq_sample, colnames(vsearchT))
setdiff(colnames(vsearchT), dukeSamples$Seq_sample)
