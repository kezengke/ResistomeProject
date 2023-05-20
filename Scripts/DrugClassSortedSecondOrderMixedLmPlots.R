#sorted second order mixed linear model plots.
rm(list = ls())
library(stringr)
library("dplyr")
library("nlme")

geneMatchT<-read.csv("MASTER_AMRlist_2023_03.csv", sep = ",", header = T)
metaData<-read.csv("metaWithBins.csv", header = T, row.names = 1)

#gene counts tables
amrT<-read.csv("CountsTables/amrFiltered.csv", header = T, row.names = 1, check.names = F)
rgiT<-read.csv("CountsTables/rgiFiltered.csv", header = T, row.names = 1, check.names = F)
vsearchT<-read.csv("CountsTables/vsearchFiltered.csv", header = T, row.names = 1, check.names = F)

# Function to find row numbers
find_row_numbers <- function(matrix, strings) {
  row_numbers <- lapply(strings, function(string) {
    matches <- which(matrix == string, arr.ind = TRUE)[, "row"]
    if (length(matches) > 0) {
      matches[1]
    } else {
      NA
    }
  })
  unlist(row_numbers)
}

amrRows<-find_row_numbers(geneMatchT, rownames(amrT))
rgiRows<-find_row_numbers(geneMatchT, rownames(rgiT))
vsearchRows<-find_row_numbers(geneMatchT, rownames(vsearchT))

#AMR
#keep only the rows present in amr counts table in master list
amrMatchT<-geneMatchT[amrRows, , drop = F]
amrMatchT<-data.frame(rownames(amrT), amrMatchT[, 12]) #keep only the gene name and Corrected AMR classification column
colnames(amrMatchT)<-c("gene", "Drug.Class")
amrMatchT<-na.omit(amrMatchT)

drugClass<-unique(amrMatchT$Drug.Class)
newamrT<-c()
for (i in 1:length(drugClass)) {
  genes<-amrMatchT %>% filter(amrMatchT$Drug.Class == drugClass[i])
  currentFrame<-amrT[genes$gene, ]
  sumRow<-colSums(currentFrame)
  newamrT<-rbind(newamrT, sumRow)
}
rownames(newamrT)<-drugClass

#RGI
#keep only the rows present in rgi counts table in master list
rgiMatchT<-geneMatchT[rgiRows, , drop = F]
rgiMatchT<-data.frame(rownames(rgiT), rgiMatchT[, 12]) #keep only the gene name and Corrected AMR classification column
colnames(rgiMatchT)<-c("gene", "Drug.Class")
rgiMatchT<-na.omit(rgiMatchT)

drugClass<-unique(rgiMatchT$Drug.Class)
newrgiT<-c()
for (i in 1:length(drugClass)) {
  genes<-rgiMatchT %>% filter(rgiMatchT$Drug.Class == drugClass[i])
  currentFrame<-rgiT[genes$gene, ]
  sumRow<-colSums(currentFrame)
  newrgiT<-rbind(newrgiT, sumRow)
}
rownames(newrgiT)<-drugClass

#VSEARCH
#keep only the rows present in vsearch counts table in master list
vsearchMatchT<-geneMatchT[vsearchRows, , drop = F]
vsearchMatchT<-data.frame(rownames(vsearchT), vsearchMatchT[, 12]) #keep only the gene name and Corrected AMR classification column
colnames(vsearchMatchT)<-c("gene", "Drug.Class")
vsearchMatchT<-na.omit(vsearchMatchT)

drugClass<-unique(vsearchMatchT$Drug.Class)
newvsearchT<-c()
for (i in 1:length(drugClass)) {
  genes<-vsearchMatchT %>% filter(vsearchMatchT$Drug.Class == drugClass[i])
  currentFrame<-vsearchT[genes$gene, ]
  sumRow<-colSums(currentFrame)
  newvsearchT<-rbind(newvsearchT, sumRow)
}
rownames(newvsearchT)<-drugClass

metaAMR<-metaData[colnames(newamrT), , drop = F]
metaRGI<-metaData[colnames(newrgiT), , drop = F]
metaVSEARCH<-metaData[colnames(newvsearchT), , drop = F]

#lm 2nd order
pdf("Plots/SortedSecondOrderMixedLm(AMRDrugClass).pdf", width=12, height=18)
par(mfrow=c(3, 2))
par(mar=c(5, 6, 4, 1)+.1)
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
pvals<-vector()
Fpval<-vector()
for (i in 1:nrow(newvsearchT)) {
  myM<-data.frame(unlist(newvsearchT[i, ]), metaVSEARCH$Timepoint, metaVSEARCH$ID)
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
  mainText<-paste0(rownames(newvsearchT)[plot_order[i]],  "\nANOVA P=", signif(pvals[plot_order[i]], 3), ", ", "adj. ANOVA-P=", signif(adjPvals[plot_order[i]], 3), 
                   "\nLine vs. Poly F-P=", signif(Fpval[plot_order[i]], 3), ", ", "adj. F-P=", signif(adjFpvals[plot_order[i]], 3))
  myM<-data.frame(unlist(newvsearchT[plot_order[i], ]), metaVSEARCH$Timepoint, metaVSEARCH$ID)
  colnames(myM)<-c("counts", "timePoint", "ID")
  Model<-lme(counts ~ poly(timePoint, 2), random = ~1 | ID, data = myM)
  plot(metaVSEARCH$Timepoint, unlist(newvsearchT[plot_order[i],]), 
       pch = 19, col = "olivedrab4", 
       main = mainText,
       xlab = "TimePoints", ylab = "Counts(Log10)")
  lines(sort(metaVSEARCH$Timepoint), fitted(Model)[order(metaVSEARCH$Timepoint)], 
        col = "grey", type = "l")
}
dev.off()
