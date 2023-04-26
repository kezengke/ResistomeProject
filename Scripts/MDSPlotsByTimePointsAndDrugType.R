#MDS plots colored by timepoints based on MASTER_AMRlist_2023_03.csv Drug Type column
rm(list = ls())
library("dplyr")
library(stringr)
library("vegan")
library("RColorBrewer")
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

#AMR
amrT<-read.delim("CountsTables/AMR_counts.tsv", sep = "\t", header = T, row.names = 2)
amrT<-amrT[,-1]
amrT<-amrT[, grepl("^D", colnames(amrT))]
colnames(amrT)<-gsub(".amrfinder.txt", "", colnames(amrT)) #get rid of useless info
SampleWithDots<-sapply(str_count(colnames(amrT), "\\."), `[`, 1) == 2
colnames(amrT)[SampleWithDots]<-sub("\\.", "", colnames(amrT)[SampleWithDots]) #get rid of random "."s in sample names
colnames(amrT)<-sub("\\.", "-", colnames(amrT)) #replace . with - to match with metadata
colnames(amrT)<-sapply(str_split(colnames(amrT), "_", n = 2), `[`, 2)
amrT<-amrT[, -27]

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

rownames(vsearchT)<-sapply(str_split(rownames(vsearchT), "\\|", n = 6), `[`, 6)

#read in gene match table
geneMatchT<-read.csv("MASTER_AMRlist_2023_03.csv", sep = ",", header = T)

#AMR
#keep only the rows present in amr counts table in master list
amrMatchT<-geneMatchT %>% filter(row_number() %in% match(rownames(amrT), geneMatchT$AMR.name))
amrMatchT<-na.omit(amrMatchT)
amrMatchT<-amrMatchT[, c(18,12)] #keep only the gene name and AMR Gene Family columns
#parse out the drug class column
drugClass<-amrMatchT$Drug.Class
drugType<-unique(unlist(str_split(drugClass, ";")))

newamrT<-c()
for (i in 1:length(drugType)) {
  allIndex<-which(sapply(str_split(drugClass, ";"), function(x) drugType[i] %in% unlist(x)))
  currentFrame<-amrT[allIndex, ]
  sumRow<-colSums(currentFrame)
  newamrT<-rbind(newamrT, sumRow)
}
rownames(newamrT)<-drugType

#RGI
#keep only the rows present in rgi counts table in master list
rgiMatchT<-geneMatchT %>% filter(row_number() %in% match(rownames(rgiT), geneMatchT$RGI.CARD.Short.Name))
rgiMatchT<-rgiMatchT[, c(15,12)] #keep only the gene name and AMR Gene Family columns

#parse out the drug class column
drugClass<-rgiMatchT$Drug.Class
drugType<-unique(unlist(str_split(drugClass, ";")))

newrgiT<-c()
for (i in 1:length(drugType)) {
  allIndex<-which(sapply(str_split(drugClass, ";"), function(x) drugType[i] %in% unlist(x)))
  currentFrame<-rgiT[allIndex, ]
  sumRow<-colSums(currentFrame)
  newrgiT<-rbind(newrgiT, sumRow)
}
rownames(newrgiT)<-drugType

#VSEARCH
#keep only the rows present in vsearch counts table in master list
vsearchMatchT<-geneMatchT %>% filter(row_number() %in% match(rownames(vsearchT), geneMatchT$Vsearch.ARO.Name))
vsearchMatchT<-vsearchMatchT[, c(7,12)] #keep only the gene name and AMR Gene Family columns

#parse out the drug class column
drugClass<-vsearchMatchT$Drug.Class
drugType<-unique(unlist(str_split(drugClass, ";")))

newvsearchT<-c()
for (i in 1:length(drugType)) {
  allIndex<-which(sapply(str_split(drugClass, ";"), function(x) drugType[i] %in% unlist(x)))
  currentFrame<-vsearchT[allIndex, ]
  sumRow<-colSums(currentFrame)
  newvsearchT<-rbind(newvsearchT, sumRow)
}
rownames(newvsearchT)<-drugType

metaAMR<-dukeSamples[colnames(newamrT), ]
metaRGI<-dukeSamples[colnames(newrgiT), ]
metaVsearch<-dukeSamples[colnames(newvsearchT), ]

pdf("Plots/MDSForDrugClass(ColoredByTimePoints).pdf", width=18, height=6)
par(mfrow=c(1,3))
par(mar=c(5,6,4,1)+.1)

circleCol<-brewer.pal(length(unique(metaAMR$bins)), "Spectral")
cols<-circleCol[factor(metaAMR$bins, levels = c("PRE", "D7", "D14", "D21", "D35", "D60", "D100"))]
MDS<-capscale(t(newamrT)~1,distance = "bray")
percentVariance<-MDS$CA$eig/sum(eigenvals(MDS))*100
pval<-adonis2(t(newamrT)~metaAMR$bins, method="bray")$aov.tab$`Pr(>F)`[1]
statusPlot<-ordiplot(MDS,choices=c(1,2),type="none",cex.lab=1,
                     xlab=paste("MDS1  ", format(percentVariance[1], digits = 4), "%", sep = ""),
                     ylab=paste("MDS2  ", format(percentVariance[2], digits = 4), "%", sep = ""),
                     main=paste("AMR(DrugType) n =", ncol(newamrT), "\nP-value:",pval))
points(statusPlot,"sites", pch=19, cex=2.5, col=adjustcolor(cols, alpha.f = 0.5))
ordiellipse(statusPlot, metaAMR$bins, kind="se", conf=0.95, lwd=4, draw = "lines", col=circleCol) 
legend("topleft", c("PRE", "D7", "D14", "D21", "D35", "D60", "D100"), col = circleCol[1:11], cex = 1, pch = 16, bty = "n")

circleCol<-brewer.pal(length(unique(metaRGI$bins)), "Spectral")
cols<-circleCol[factor(metaRGI$bins, levels = c("PRE", "D7", "D14", "D21", "D35", "D60", "D100"))]
MDS<-capscale(t(newrgiT)~1,distance = "bray")
percentVariance<-MDS$CA$eig/sum(eigenvals(MDS))*100
pval<-adonis2(t(newrgiT)~metaRGI$bins, method="bray")$aov.tab$`Pr(>F)`[1]
statusPlot<-ordiplot(MDS,choices=c(1,2),type="none",cex.lab=1,
                     xlab=paste("MDS1  ", format(percentVariance[1], digits = 4), "%", sep = ""),
                     ylab=paste("MDS2  ", format(percentVariance[2], digits = 4), "%", sep = ""),
                     main=paste("RGI(DrugType) n =", ncol(newrgiT), "\nP-value:",pval))
points(statusPlot,"sites", pch=19, cex=2.5, col=adjustcolor(cols, alpha.f = 0.5))
ordiellipse(statusPlot, metaRGI$bins, kind="se", conf=0.95, lwd=4, draw = "lines", col=circleCol) 
legend("topleft", c("PRE", "D7", "D14", "D21", "D35", "D60", "D100"), col = circleCol[1:11], cex = 1, pch = 16, bty = "n")

circleCol<-brewer.pal(length(unique(metaVsearch$bins)), "Spectral")
cols<-circleCol[factor(metaVsearch$bins, levels = c("PRE", "D7", "D14", "D21", "D35", "D60", "D100"))]
MDS<-capscale(t(newvsearchT)~1,distance = "bray")
percentVariance<-MDS$CA$eig/sum(eigenvals(MDS))*100
pval<-adonis2(t(newvsearchT)~metaVsearch$bins, method="bray")$aov.tab$`Pr(>F)`[1]
statusPlot<-ordiplot(MDS,choices=c(1,2),type="none",cex.lab=1,
                     xlab=paste("MDS1  ", format(percentVariance[1], digits = 4), "%", sep = ""),
                     ylab=paste("MDS2  ", format(percentVariance[2], digits = 4), "%", sep = ""),
                     main=paste("vsearch(DrugType) n =", ncol(newvsearchT), "\nP-value:",pval))
points(statusPlot,"sites", pch=19, cex=2.5, col=adjustcolor(cols, alpha.f = 0.5))
ordiellipse(statusPlot, metaVsearch$bins, kind="se", conf=0.95, lwd=4, draw = "lines", col=circleCol) 
legend("topright", c("PRE", "D7", "D14", "D21", "D35", "D60", "D100"), col = circleCol[1:11], cex = 1, pch = 16, bty = "n")

dev.off()