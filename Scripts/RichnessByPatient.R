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

# #Normalization
# n<-colSums(brackenT)
# sumx<-sum(brackenT)
# brackenT<-log10((brackenT/n)*(sumx/ncol(brackenT))+1)

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

pdf("Plots/TypeRichnessByPatients.pdf", width=12, height=12)
par(mfrow=c(2, 2))
par(mar=c(5, 6, 4, 1)+.1)

binCounts<-aggregate(colSums(brackenT>0), list(dukeSamples[colnames(brackenT), 5]), FUN=sum)
barplot(binCounts$x, names.arg=binCounts$Group.1, las = 3,
        main = "Taxa Richness by Patient (Bracken-Species)", col = adjustcolor("tan2", alpha.f = 0.7),
        ylab = "Type of Taxa",
        cex.axis = 1, cex.lab = 1.5, cex.names = 1, cex.main = 1.8)

binCounts<-aggregate(colSums(amrT>0), list(dukeSamples[colnames(amrT), 5]), FUN=sum)
barplot(binCounts$x, names.arg=binCounts$Group.1, las = 3, 
        main = "Gene Richness by Patient (AMR)", col = adjustcolor("coral3", alpha.f = 0.7), 
        ylab = "Type of Genes",
        cex.axis = 1, cex.lab = 1.5, cex.names = 1, cex.main = 1.8)

binCounts<-aggregate(colSums(rgiT>0), list(dukeSamples[colnames(rgiT), 5]), FUN=sum)
barplot(binCounts$x, names.arg=binCounts$Group.1, las = 3, 
        main = "Gene Richness by Patient (RGI)", col = adjustcolor("cornflowerblue", alpha.f = 0.7),
        ylab = "Type of Genes", 
        cex.axis = 1, cex.lab = 1.5, cex.names = 1, cex.main = 1.8)

binCounts<-aggregate(colSums(vsearchT>0), list(dukeSamples[colnames(vsearchT), 5]), FUN=sum)
barplot(binCounts$x, names.arg=binCounts$Group.1, las = 3,
        main = "Gene Richness by Patient (vsearch)", col = adjustcolor("olivedrab4", alpha.f = 0.7),
        ylab = "Type of Genes", 
        cex.axis = 1, cex.lab = 1.5, cex.names = 1, cex.main = 1.8,
)

dev.off()
