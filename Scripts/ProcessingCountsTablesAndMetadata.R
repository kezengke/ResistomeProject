#process metadata, bracken and three gene counts tables.
rm(list = ls())
library(stringr)
library("dplyr")

#proposed bins
proposedBins<-read.csv("ProposedBin.csv", header = T, row.names = 1)
rownames(proposedBins)<-sapply(str_split(rownames(proposedBins), "_", n = 2), `[`, 2)
proposedBins$Proposed.bin[proposedBins$Proposed.bin != "PRE"]<-ifelse(as.numeric(proposedBins$Proposed.bin[proposedBins$Proposed.bin != "PRE"]) >= 180, 
                                  180, 
                                  proposedBins$Proposed.bin[proposedBins$Proposed.bin != "PRE"])
proposedBins$Proposed.bin<-ifelse(substr(proposedBins$Proposed.bin, 1, 3) == "PRE",
                                  proposedBins$Proposed.bin,
                                  paste0("D", proposedBins$Proposed.bin))

#in and out patient
ptInOut<-read.csv("patientInOut.csv", header = T, row.names = 1)

#metadata
dukeSamples<-read.csv("Duke_samples_meta.csv", header = T, sep = ",")
rownames(dukeSamples)<-dukeSamples$Seq_sample
rownames(dukeSamples)<-sapply(str_split(rownames(dukeSamples), "_", n = 2), `[`, 2)
dukeSamples$ID<-paste0("D", gsub("pre","" , sapply(str_split(dukeSamples$sample, "D", n = 3), `[`, 2), ignore.case = T)) #adding sample ID

dukeSamples$bins<-proposedBins$Proposed.bin

dukeSamples$ptInOut<-ptInOut[as.character(dukeSamples$PID),]

write.csv(dukeSamples, "metaWithBins.csv")

#bracken
brackenT<-read.delim("CountsTables/duke_bracken.csv", sep = ",", header = T, row.names = 1)
brackenT<-brackenT[, -c(1,2)] #get rid of taxonomy ID and level
brackenT<-brackenT[, grepl("num", colnames(brackenT))] #get rid of fractions

colnames(brackenT)<-gsub(".bracken.out_num", "", colnames(brackenT)) #get rid of useless info
SampleWithDots<-sapply(str_count(colnames(brackenT), "\\."), `[`, 1) == 2
colnames(brackenT)[SampleWithDots]<-sub("\\.", "", colnames(brackenT)[SampleWithDots]) #get rid of random "."s in sample names
colnames(brackenT)<-sub("\\.", "-", colnames(brackenT)) #replace . with - to match with metadata
colnames(brackenT)<-sapply(str_split(colnames(brackenT), "_", n = 2), `[`, 2) #keeping sequencing info for matching

write.csv(brackenT, "CountsTables/brackenProcessed.csv")

#AMR
amrT<-read.delim("CountsTables/AMR_counts.tsv", sep = "\t", header = T, row.names = 2)
amrT<-amrT[,-1]
amrT<-amrT[, grepl("^D", colnames(amrT))]
colnames(amrT)<-gsub(".amrfinder.txt", "", colnames(amrT)) #get rid of useless info
SampleWithDots<-sapply(str_count(colnames(amrT), "\\."), `[`, 1) == 2
colnames(amrT)[SampleWithDots]<-sub("\\.", "", colnames(amrT)[SampleWithDots]) #get rid of random "."s in sample names
colnames(amrT)<-sub("\\.", "-", colnames(amrT)) #replace . with - to match with metadata
colnames(amrT)<-sapply(str_split(colnames(amrT), "_", n = 2), `[`, 2)
#get rid of sample with all 0s
amrT<-amrT[, -27]

write.csv(amrT, "CountsTables/amrProcessed.csv")

#RGI
rgiT<-read.delim("CountsTables/RGI_counts.tsv", sep = "\t", header = T, row.names = 2)
rgiT<-rgiT[, -1]
rgiT<-rgiT[, grepl("^D", colnames(rgiT))]
colnames(rgiT)<-gsub(".rgi.txt", "", colnames(rgiT)) #get rid of useless info
SampleWithDots<-sapply(str_count(colnames(rgiT), "\\."), `[`, 1) == 2
colnames(rgiT)[SampleWithDots]<-sub("\\.", "", colnames(rgiT)[SampleWithDots]) #get rid of random "."s in sample names
colnames(rgiT)<-sub("\\.", "-", colnames(rgiT)) #replace . with - to match with metadata
colnames(rgiT)<-sapply(str_split(colnames(rgiT), "_", n = 2), `[`, 2)

write.csv(rgiT, "CountsTables/rgiProcessed.csv")

#vsearch
vsearchT<-read.delim("CountsTables/vsearch_counts.tsv", sep = "\t", header = T, row.names = 2)
vsearchT<-vsearchT[, -1] #get rid of index column
vsearchT<-vsearchT[, grepl("^D", colnames(vsearchT))]
colnames(vsearchT)<-gsub(".txt", "", colnames(vsearchT)) #get rid of useless inåfo

SampleWithDots<-sapply(str_count(colnames(vsearchT), "\\."), `[`, 1) == 2
colnames(vsearchT)[SampleWithDots]<-sub("\\.", "", colnames(vsearchT)[SampleWithDots]) #get rid of random "."s in sample names
colnames(vsearchT)<-sub("\\.", "-", colnames(vsearchT)) #replace . with - to match with metadata
colnames(vsearchT)<-sapply(str_split(colnames(vsearchT), "_", n = 2), `[`, 2)

rownames(vsearchT)<-sapply(str_split(rownames(vsearchT), "\\|", n = 6), `[`, 6)

write.csv(vsearchT, "CountsTables/vsearchProcessed.csv")