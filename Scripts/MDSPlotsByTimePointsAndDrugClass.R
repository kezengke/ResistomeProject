#MDS plots colored by timepoints based on MASTER_AMRlist_2023_03.csv Drug Type column
rm(list = ls())
library("vegan")
library(stringr)
library("dplyr")
library("RColorBrewer")

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

pdf("Plots/MDSForDrugClass(ColoredByTimePoints).pdf", width=18, height=6)
par(mfrow=c(1,3))
par(mar=c(5,6,4,1)+.1)

circleCol<-brewer.pal(length(unique(metaAMR$bins)), "Spectral")
cols<-circleCol[factor(metaAMR$bins, levels = c("PRE", "D0", "D7", "D14", "D21", "D28", "D35", "D60", "D100", "D180"))]
MDS<-capscale(t(newamrT)~1,distance = "bray")
percentVariance<-MDS$CA$eig/sum(eigenvals(MDS))*100
pval<-adonis2(t(newamrT)~metaAMR$bins, method="bray")$`Pr(>F)`[1]
statusPlot<-ordiplot(MDS,choices=c(1,2),type="none",cex.lab=1,
                     xlab=paste("MDS1  ", format(percentVariance[1], digits = 4), "%", sep = ""),
                     ylab=paste("MDS2  ", format(percentVariance[2], digits = 4), "%", sep = ""),
                     main=paste("AMR(DrugClass) n =", ncol(newamrT), "\nP-value:",pval))
points(statusPlot,"sites", pch=19, cex=2.5, col=adjustcolor(cols, alpha.f = 0.5))
for (i in 1:length(unique(metaAMR$bins))) {
  ordiellipse(statusPlot, metaAMR$bins, kind="se", conf=0.95, lwd=4, draw = "lines", col=circleCol[i], 
              show.groups = c("PRE", "D0", "D7", "D14", "D21", "D28", "D35", "D60", "D100", "D180")[i],
              font=2,cex=1) 
}
legend("topleft", c("PRE", "D0", "D7", "D14", "D21", "D28", "D35", "D60", "D100", "D180"), col = circleCol[1:10], cex = 1, pch = 16, bty = "n")

circleCol<-brewer.pal(length(unique(metaRGI$bins)), "Spectral")
cols<-circleCol[factor(metaRGI$bins, levels = c("PRE", "D0", "D7", "D14", "D21", "D28", "D35", "D60", "D100", "D180"))]
MDS<-capscale(t(newrgiT)~1,distance = "bray")
percentVariance<-MDS$CA$eig/sum(eigenvals(MDS))*100
pval<-adonis2(t(newrgiT)~metaRGI$bins, method="bray")$`Pr(>F)`[1]
statusPlot<-ordiplot(MDS,choices=c(1,2),type="none",cex.lab=1,
                     xlab=paste("MDS1  ", format(percentVariance[1], digits = 4), "%", sep = ""),
                     ylab=paste("MDS2  ", format(percentVariance[2], digits = 4), "%", sep = ""),
                     main=paste("RGI(DrugClass) n =", ncol(newrgiT), "\nP-value:",pval))
points(statusPlot,"sites", pch=19, cex=2.5, col=adjustcolor(cols, alpha.f = 0.5))
for (i in 1:length(unique(metaRGI$bins))) {
  ordiellipse(statusPlot, metaRGI$bins, kind="se", conf=0.95, lwd=4, draw = "lines", col=circleCol[i], 
              show.groups = c("PRE", "D0", "D7", "D14", "D21", "D28", "D35", "D60", "D100", "D180")[i],
              font=2,cex=1) 
}
legend("topleft", c("PRE", "D0", "D7", "D14", "D21", "D28", "D35", "D60", "D100", "D180"), col = circleCol[1:10], cex = 1, pch = 16, bty = "n")

circleCol<-brewer.pal(length(unique(metaVSEARCH$bins)), "Spectral")
cols<-circleCol[factor(metaVSEARCH$bins, levels = c("PRE", "D0", "D7", "D14", "D21", "D28", "D35", "D60", "D100", "D180"))]
MDS<-capscale(t(newvsearchT)~1,distance = "bray")
percentVariance<-MDS$CA$eig/sum(eigenvals(MDS))*100
pval<-adonis2(t(newvsearchT)~metaVSEARCH$bins, method="bray")$`Pr(>F)`[1]
statusPlot<-ordiplot(MDS,choices=c(1,2),type="none",cex.lab=1,
                     xlab=paste("MDS1  ", format(percentVariance[1], digits = 4), "%", sep = ""),
                     ylab=paste("MDS2  ", format(percentVariance[2], digits = 4), "%", sep = ""),
                     main=paste("vsearch(DrugClass) n =", ncol(newvsearchT), "\nP-value:",pval))
points(statusPlot,"sites", pch=19, cex=2.5, col=adjustcolor(cols, alpha.f = 0.5))
for (i in 1:length(unique(metaVSEARCH$bins))) {
  ordiellipse(statusPlot, metaVSEARCH$bins, kind="se", conf=0.95, lwd=4, draw = "lines", col=circleCol[i], 
              show.groups = c("PRE", "D0", "D7", "D14", "D21", "D28", "D35", "D60", "D100", "D180")[i],
              font=2,cex=1) 
} 
legend("topright", c("PRE", "D0", "D7", "D14", "D21", "D28", "D35", "D60", "D100", "D180"), col = circleCol[1:10], cex = 1, pch = 16, bty = "n")

dev.off()