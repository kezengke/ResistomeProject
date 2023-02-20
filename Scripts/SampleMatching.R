rm(list = ls())
library(stringr)
library("dplyr")
dukeSamples<-read.csv("Duke_samples_meta.csv", header = T, sep = ",")
rownames(dukeSamples)<-dukeSamples$Seq_sample
rownames(dukeSamples)<-sapply(str_split(rownames(dukeSamples), "_", n = 2), `[`, 2)

matchingT<-dukeSamples
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
brackenTag<-colnames(brackenT)
brackenTag<-gsub(".bracken.out_num", "", brackenTag) #get rid of useless info
SampleWithDots<-sapply(str_count(brackenTag, "\\."), `[`, 1) == 2
brackenTag[SampleWithDots]<-sub("\\.", "", brackenTag[SampleWithDots]) #get rid of random "."s in sample names
brackenTag<-sub("\\.", "-", brackenTag) #replace . with - to match with metadata
brackenTag<-sapply(str_split(brackenTag, "_", n = 2), `[`, 2) #keeping sequencing info for matching

brackenNames<-data.frame(colnames(brackenT))
rownames(brackenNames)<-brackenTag

matchingT$brackenSampleNames<-brackenNames[rownames(matchingT), ]

#AMR
amrT<-read.delim("CountsTables/AMR_counts.tsv", sep = "\t", header = T, row.names = 1)
amrT<-amrT[, grepl("^D", colnames(amrT))]
amrTag<-colnames(amrT)
amrTag<-gsub(".amrfinder.txt", "", amrTag) #get rid of useless info
SampleWithDots<-sapply(str_count(amrTag, "\\."), `[`, 1) == 2
amrTag[SampleWithDots]<-sub("\\.", "", amrTag[SampleWithDots]) #get rid of random "."s in sample names
amrTag<-sub("\\.", "-", amrTag) #replace . with - to match with metadata
amrTag<-sapply(str_split(amrTag, "_", n = 2), `[`, 2) #get rid of the potentially mislabeled portion

amrNames<-data.frame(colnames(amrT))
rownames(amrNames)<-amrTag

matchingT$AMRsampleNames<-amrNames[rownames(matchingT), ]

#RGI
rgiT<-read.delim("CountsTables/RGI_counts.tsv", sep = "\t", header = T, row.names = 1)
rgiT<-rgiT[, grepl("^D", colnames(rgiT))]
rgiTag<-colnames(rgiT)
rgiTag<-gsub(".rgi.txt", "", rgiTag) #get rid of useless info
SampleWithDots<-sapply(str_count(rgiTag, "\\."), `[`, 1) == 2
rgiTag[SampleWithDots]<-sub("\\.", "", rgiTag[SampleWithDots]) #get rid of random "."s in sample names
rgiTag<-sub("\\.", "-", rgiTag) #replace . with - to match with metadata
rgiTag<-sapply(str_split(rgiTag, "_", n = 2), `[`, 2)

rgiNames<-data.frame(colnames(rgiT))
rownames(rgiNames)<-rgiTag

matchingT$RGIsampleNames<-rgiNames[rownames(matchingT), ]

#vsearch
vsearchT<-read.delim("CountsTables/vsearch_counts.tsv", sep = "\t", header = T, row.names = 2)
vsearchT<-vsearchT[, -1] #get rid of index column
vsearchT<-vsearchT[, grepl("^D", colnames(vsearchT))]
vsearchTag<-colnames(vsearchT)
vsearchTag<-gsub(".txt", "", vsearchTag) #get rid of useless info
SampleWithDots<-sapply(str_count(vsearchTag, "\\."), `[`, 1) == 2
vsearchTag[SampleWithDots]<-sub("\\.", "", vsearchTag[SampleWithDots]) #get rid of random "."s in sample names
vsearchTag<-sub("\\.", "-", vsearchTag) #replace . with - to match with metadata
vsearchTag<-sapply(str_split(vsearchTag, "_", n = 2), `[`, 2)

vsearchNames<-data.frame(colnames(vsearchT))
rownames(vsearchNames)<-vsearchTag

matchingT$VsearchSampleNames<-vsearchNames[rownames(matchingT), ]


write.table(matchingT, "SampleMatchingTable.txt", sep = "\t", row.names = F)
