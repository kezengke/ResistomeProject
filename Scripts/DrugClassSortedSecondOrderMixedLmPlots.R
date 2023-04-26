rm(list = ls())
library(stringr)
library("RColorBrewer")
library("vegan")
library("dplyr")
library("nlme")

#gene category match table
geneMatchT<-read.csv("MASTER_AMRlist_2023_03.csv", sep = ",", header = T)
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

#AMR
#keep only the rows present in amr counts table in master list
amrMatchT<-geneMatchT %>% filter(row_number() %in% match(rownames(amrT), geneMatchT$AMR.name))
amrMatchT<-na.omit(amrMatchT)
amrMatchT<-amrMatchT[, c(18,12)] #keep only the gene name and AMR classification columns

amrTypes<-unique(amrMatchT$Drug.Class) #all drug types
newamrT<-c()
for (i in 1:length(amrTypes)) {
  genes<-amrMatchT %>% filter(amrMatchT$Drug.Class == amrTypes[i])
  currentFrame<-amrT[genes$AMR.name, ]
  sumRow<-colSums(currentFrame)
  newamrT<-rbind(newamrT, sumRow)
}
rownames(newamrT)<-amrTypes
newamrT<-newamrT[, -which(colSums(newamrT) == 0)]

#RGI
#keep only the rows present in rgi counts table in master list
rgiMatchT<-geneMatchT %>% filter(row_number() %in% match(rownames(rgiT), geneMatchT$RGI.CARD.Short.Name))
rgiMatchT<-rgiMatchT[, c(15,12)] #keep only the gene name and AMR classification columns

rgiTypes<-unique(rgiMatchT$Drug.Class) #all drug types
newrgiT<-c()
for (i in 1:length(rgiTypes)) {
  genes<-rgiMatchT %>% filter(rgiMatchT$Drug.Class == rgiTypes[i])
  currentFrame<-rgiT[genes$RGI.CARD.Short.Name, ]
  sumRow<-colSums(currentFrame)
  newrgiT<-rbind(newrgiT, sumRow)
}
rownames(newrgiT)<-rgiTypes

#VSEARCH
#keep only the rows present in vsearch counts table in master list
vsearchMatchT<-geneMatchT %>% filter(row_number() %in% match(rownames(vsearchT), geneMatchT$Vsearch.ARO.Name))
vsearchMatchT<-vsearchMatchT[, c(7,12)] #keep only the gene name and AMR classification columns

vsearchTypes<-unique(vsearchMatchT$Drug.Class) #all drug types
newvsearchT<-c()
for (i in 1:length(vsearchTypes)) {
  genes<-vsearchMatchT %>% filter(vsearchMatchT$Drug.Class == vsearchTypes[i])
  currentFrame<-vsearchT[genes$Vsearch.ARO.Name, ]
  sumRow<-colSums(currentFrame)
  newvsearchT<-rbind(newvsearchT, sumRow)
}
rownames(newvsearchT)<-vsearchTypes

#meta for each table
metaAMR<-dukeSamples[colnames(newamrT), ]
metaRGI<-dukeSamples[colnames(newrgiT), ]
metaVsearch<-dukeSamples[colnames(newvsearchT), ]

#Normalization
n<-colSums(newamrT)
sumx<-sum(newamrT)
for (i in 1:ncol(newamrT)) {
  newamrT[,i]<-newamrT[,i]/n[i]
}
newamrT<-log10(newamrT*(sumx/ncol(newamrT))+1)

n<-colSums(newrgiT)
sumx<-sum(newrgiT)
for (i in 1:ncol(newrgiT)) {
  newrgiT[,i]<-newrgiT[,i]/n[i]
}
newrgiT<-log10(newrgiT*(sumx/ncol(newrgiT))+1)

n<-colSums(newvsearchT)
sumx<-sum(newvsearchT)
for (i in 1:ncol(newvsearchT)) {
  newvsearchT[,i]<-newvsearchT[,i]/n[i]
}
newvsearchT<-log10(newvsearchT*(sumx/ncol(newvsearchT))+1)

#lm 2nd order
pdf("Plots/SortedSecondOrderMixedLm(AMRDrugClass).pdf", width=12, height=18)
par(mfrow=c(3, 2))
par(mar=c(5, 6, 4, 1)+.1)
newamrT<-newamrT[apply(newamrT == 0, 1, sum) <= (ncol(newamrT)*0.8), ]
pvals<-vector()
Fpval<-vector()
for (i in 1:nrow(newamrT)) {
  myM<-data.frame(unlist(newamrT[i, ]), metaAMR$Timepoint, metaAMR$ID)
  colnames(myM)<-c("counts", "timePoint", "ID")
  lineM<-lme(counts ~ timePoint, random = ~1 | ID, data = myM)
  polyM<-lme(counts ~ poly(timePoint, 2), random = ~1 | ID, data = myM)
  
  reducedError<-sum(resid(lineM)^2)
  fullError<-sum(resid(polyM)^2)
  
  reducedDF<-nrow(myM)-2
  fullDF<-nrow(myM)-3
  
  myF <- ((reducedError - fullError)/(reducedDF - fullDF))/(fullError/fullDF)
  
  Fpval[i]<-pf(myF, 1, fullDF, lower.tail = F)
  pvals[i]<-anova(polyM)[2,4]
}
adjPvals<-p.adjust(pvals, method = "BH")
adjFpvals<-p.adjust(Fpval, method = "BH")

plot_order<-order(pvals)
for (i in 1:nrow(newamrT)) {
  mainText<-paste0(rownames(newamrT)[plot_order[i]], "\nANOVA P=", signif(pvals[plot_order[i]], 3), ", ", "adj. ANOVA-P=", signif(adjPvals[plot_order[i]], 3), 
                   "\nLine vs. Poly F-P=", signif(Fpval[plot_order[i]], 3), ", ", "adj. F-P=", signif(adjFpvals[plot_order[i]], 3))
  myM<-data.frame(unlist(newamrT[plot_order[i], ]), metaAMR$Timepoint, metaAMR$ID)
  colnames(myM)<-c("counts", "timePoint", "ID")
  Model<-lme(counts ~ poly(timePoint, 2), random = ~1 | ID, data = myM)
  plot(metaAMR$Timepoint, unlist(newamrT[plot_order[i],]), 
       pch = 19, col = "coral3", 
       main = mainText,
       xlab = "TimePoints", ylab = "Counts(Log10)")
  lines(sort(metaAMR$Timepoint), fitted(Model)[order(metaAMR$Timepoint)], 
        col = "grey", type = "l")
}
dev.off()

pdf("Plots/SortedSecondOrderMixedLm(RGIDrugClass).pdf", width=12, height=18)
par(mfrow=c(3, 2))
par(mar=c(5, 6, 4, 1)+.1)
newrgiT<-newrgiT[apply(newrgiT == 0, 1, sum) <= (ncol(newrgiT)*0.8), ]
pvals<-vector()
Fpval<-vector()
for (i in 1:nrow(newrgiT)) {
  myM<-data.frame(unlist(newrgiT[i, ]), metaRGI$Timepoint, metaRGI$ID)
  colnames(myM)<-c("counts", "timePoint", "ID")
  lineM<-lme(counts ~ timePoint, random = ~1 | ID, data = myM)
  polyM<-lme(counts ~ poly(timePoint, 2), random = ~1 | ID, data = myM)
  
  reducedError<-sum(resid(lineM)^2)
  fullError<-sum(resid(polyM)^2)
  
  reducedDF<-nrow(myM)-2
  fullDF<-nrow(myM)-3
  
  myF <- ((reducedError - fullError)/(reducedDF - fullDF))/(fullError/fullDF)
  
  Fpval[i]<-pf(myF, 1, fullDF, lower.tail = F)
  pvals[i]<-anova(polyM)[2,4]
}
adjPvals<-p.adjust(pvals, method = "BH")
adjFpvals<-p.adjust(Fpval, method = "BH")

plot_order<-order(pvals)
for (i in 1:nrow(newrgiT)) {
  mainText<-paste0(rownames(newrgiT)[plot_order[i]], "\nANOVA P=", signif(pvals[plot_order[i]], 3), ", ", "adj. ANOVA-P=", signif(adjPvals[plot_order[i]], 3), 
                   "\nLine vs. Poly F-P=", signif(Fpval[plot_order[i]], 3), ", ", "adj. F-P=", signif(adjFpvals[plot_order[i]], 3))
  myM<-data.frame(unlist(newrgiT[plot_order[i], ]), metaRGI$Timepoint, metaRGI$ID)
  colnames(myM)<-c("counts", "timePoint", "ID")
  Model<-lme(counts ~ poly(timePoint, 2), random = ~1 | ID, data = myM)
  plot(metaRGI$Timepoint, unlist(newrgiT[plot_order[i],]), 
       pch = 19, col = "cornflowerblue", 
       main = mainText,
       xlab = "TimePoints", ylab = "Counts(Log10)")
  lines(sort(metaRGI$Timepoint), fitted(Model)[order(metaRGI$Timepoint)], 
        col = "grey", type = "l")
}
dev.off()

pdf("Plots/SortedSecondOrderMixedLm(vsearchDrugClass).pdf", width=12, height=18)
par(mfrow=c(3, 2))
par(mar=c(5, 6, 4, 1)+.1)
newvsearchT<-newvsearchT[apply(newvsearchT == 0, 1, sum) <= (ncol(newvsearchT)*0.8), ]
pvals<-vector()
Fpval<-vector()
for (i in 1:nrow(newvsearchT)) {
  myM<-data.frame(unlist(newvsearchT[i, ]), metaVsearch$Timepoint, metaVsearch$ID)
  colnames(myM)<-c("counts", "timePoint", "ID")
  lineM<-lme(counts ~ timePoint, random = ~1 | ID, data = myM)
  polyM<-lme(counts ~ poly(timePoint, 2), random = ~1 | ID, data = myM)
  
  reducedError<-sum(resid(lineM)^2)
  fullError<-sum(resid(polyM)^2)
  
  reducedDF<-nrow(myM)-2
  fullDF<-nrow(myM)-3
  
  myF <- ((reducedError - fullError)/(reducedDF - fullDF))/(fullError/fullDF)
  
  Fpval[i]<-pf(myF, 1, fullDF, lower.tail = F)
  pvals[i]<-anova(polyM)[2,4]
}
adjPvals<-p.adjust(pvals, method = "BH")
adjFpvals<-p.adjust(Fpval, method = "BH")

plot_order<-order(pvals)
for (i in 1:nrow(newvsearchT)) {
  i=1
  mainText<-paste0(rownames(newvsearchT)[plot_order[i]],  "\nANOVA P=", signif(pvals[plot_order[i]], 3), ", ", "adj. ANOVA-P=", signif(adjPvals[plot_order[i]], 3), 
                   "\nLine vs. Poly F-P=", signif(Fpval[plot_order[i]], 3), ", ", "adj. F-P=", signif(adjFpvals[plot_order[i]], 3))
  myM<-data.frame(unlist(newvsearchT[plot_order[i], ]), metaVsearch$Timepoint, metaVsearch$ID)
  colnames(myM)<-c("counts", "timePoint", "ID")
  Model<-lme(counts ~ poly(timePoint, 2), random = ~1 | ID, data = myM)
  plot(metaVsearch$Timepoint, unlist(newvsearchT[plot_order[i],]), 
       pch = 19, col = "olivedrab4", 
       main = mainText,
       xlab = "TimePoints", ylab = "Counts(Log10)")
  lines(sort(metaVsearch$Timepoint), fitted(Model)[order(metaVsearch$Timepoint)], 
        col = "grey", type = "l")
}
dev.off()
