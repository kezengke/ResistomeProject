# pval histogram of mixed linear models
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

pdf("Plots/MixedSecondOrderLmDrugClassPvalHistograms.pdf", width=18, height=6)
par(mfrow=c(1, 3))
par(mar=c(5, 6, 4, 1)+.1)

modelPvals<-vector()
for (i in 1:nrow(newamrT)) {
  myM<-data.frame(unlist(newamrT[i, ]), metaAMR$Timepoint, metaAMR$ID)
  colnames(myM)<-c("counts", "timePoint", "ID")
  Model<-lme(counts ~ poly(timePoint, 2), random = ~1 | ID, data = myM)
  modelPvals[i]<-anova(Model)[2,4]
  
}
hist(modelPvals, breaks=seq(0, 1, 0.05), xlab = "p-value", col = "coral3",
     main = "AMR(DrugClass) Mixed LM 2nd Order ANOVA P-vals", 
     xlim=c(0,1), cex.lab = 1.5, cex.main = 1.7, cex.axis = 1.4)
abline(v=0.05, col=gray(.5), lty=2)
mtext(at=0.05, side=3, text=0.05, col=gray(.5))

modelPvals<-vector()
for (i in 1:nrow(newrgiT)) {
  myM<-data.frame(unlist(newrgiT[i, ]), metaRGI$Timepoint, metaRGI$ID)
  colnames(myM)<-c("counts", "timePoint", "ID")
  Model<-lme(counts ~ poly(timePoint, 2), random = ~1 | ID, data = myM)
  modelPvals[i]<-anova(Model)[2,4]
  
}
hist(modelPvals, breaks=seq(0, 1, 0.05), xlab = "p-value", col = "cornflowerblue",
     main = "RGI(DrugClass) Mixed LM 2nd Order ANOVA P-vals", 
     xlim=c(0,1), cex.lab = 1.5, cex.main = 1.7, cex.axis = 1.4)
abline(v=0.05, col=gray(.5), lty=2)
mtext(at=0.05, side=3, text=0.05, col=gray(.5))

modelPvals<-vector()
for (i in 1:nrow(newvsearchT)) {
  myM<-data.frame(unlist(newvsearchT[i, ]), metaVSEARCH$Timepoint, metaVSEARCH$ID)
  colnames(myM)<-c("counts", "timePoint", "ID")
  Model<-lme(counts ~ poly(timePoint, 2), random = ~1 | ID, data = myM)
  modelPvals[i]<-anova(Model)[2,4]
  
}
hist(modelPvals, breaks=seq(0, 1, 0.05), xlab = "p-value", col = "olivedrab4",
     main = "vsearch(DrugClass) Mixed LM 2nd Order ANOVA P-vals", 
     xlim=c(0,1), cex.lab = 1.5, cex.main = 1.7, cex.axis = 1.4)
abline(v=0.05, col=gray(.5), lty=2)
mtext(at=0.05, side=3, text=0.05, col=gray(.5))

dev.off()
