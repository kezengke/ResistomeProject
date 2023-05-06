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
for (i in 1:ncol(brackenT)) {
  brackenT[,i]<-brackenT[,i]/n[i]
}
brackenT<-log10(brackenT*(sumx/ncol(brackenT))+1)
# brackenT<-data.frame(brackenT)

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

pdf("Plots/MDS1BoxPlot(ByPatients).pdf", width=12, height=12)
par(mfrow=c(2,2))
par(mar=c(5,6,4,1)+.1)

#Bracken
MDS<-summary(capscale(t(brackenT)~1,distance = "bray"))[["sites"]][,c(1,2)]
pval<-summary(aov(MDS[,1]~metaBracken$ID))[[1]][["Pr(>F)"]][1]
boxplot(MDS~metaBracken$ID, outline = F, col = "tan2", las = 3,
        ylab = "MDS1", xlab = "", main = paste0("Bracken-Species \nANOVA pvalue:", signif(pval, digits = 3)))
# stripchart(MDS~metaBracken$ID, method = "jitter", vertical = T, pch = 19, add = T)

#AMR
MDS<-summary(capscale(t(amrT)~1,distance = "bray"))[["sites"]][,c(1,2)]
pval<-summary(aov(MDS[,1]~metaAMR$ID))[[1]][["Pr(>F)"]][1]
boxplot(MDS~metaAMR$ID, outline = F, col = "coral3", las = 3,
        ylab = "MDS1", xlab = "", main = paste0("AMR \nANOVA pvalue:", signif(pval, digits = 3)))

#RGI
MDS<-summary(capscale(t(rgiT)~1,distance = "bray"))[["sites"]][,c(1,2)]
pval<-summary(aov(MDS[,1]~metaRGI$ID))[[1]][["Pr(>F)"]][1]
boxplot(MDS~metaRGI$ID, outline = F, col = "cornflowerblue", las = 3,
        ylab = "MDS1", xlab = "", main = paste0("RGI \nANOVA pvalue:", signif(pval, digits = 3)))

#vsearch
MDS<-summary(capscale(t(vsearchT)~1,distance = "bray"))[["sites"]][,c(1,2)]
pval<-summary(aov(MDS[,1]~metaVsearch$ID))[[1]][["Pr(>F)"]][1]
boxplot(MDS~metaVsearch$ID, outline = F, col = "olivedrab4", las = 3,
        ylab = "MDS1", xlab = "", main = paste0("vsearch \nANOVA pvalue:",signif(pval, digits = 3)))

dev.off()
