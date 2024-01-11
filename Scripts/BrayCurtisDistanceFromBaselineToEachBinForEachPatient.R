#bray-curtis distance from pre to other bin for each patient
rm(list = ls())
library("vegan")
library(dplyr)
library(ggplot2)
library(stringr)

metaData<-read.csv("MetaBinnedWithInOutDonorTreatmentType.csv", header = T, row.names = 1)
#counts tables
brackenT<-read.csv("CountsTables/brackenProcessed.csv", header = T, row.names = 1, check.names = F)
brackenT<-brackenT[, intersect(rownames(metaData), colnames(brackenT)), drop = F]

distanceMatrix<-function(countsT, dissimilarityMetric){
  distanceT = as.matrix(vegdist(t(countsT), method = dissimilarityMetric))
  return(distanceT)
}

natural_sort <- function(x) {
  numeric_part <- as.numeric(gsub("\\D", "", x))
  ordered_indices <- order(numeric_part)
  x[ordered_indices]
}

brackenM<-distanceMatrix(brackenT, "bray")
names<-paste0(metaData[colnames(brackenM), 5], "_",metaData[colnames(brackenM), 6])
colnames(brackenM)<-names
rownames(brackenM)<-names

patient_ids<-unique(metaData$ID)
bins<-unique(metaData$bins[metaData$bins != "PRE"])

pdf("Plots/BrayCurtisDistanceFromBaselineForEachPatient.pdf", width=30, height=15)
par(mfrow=c(3, 5))
par(mar=c(5, 6, 4, 1)+.1)

plot(1, type = "n", xlim = c(1, length(bins)), ylim = range(brackenM),
     xlab = "Date Bin", ylab = "Bray-Curtis Distance", axes = FALSE)
axis(1, at = 1:length(bins), labels = natural_sort(bins))
axis(2)  # Add y-axis
box()
legend("topright", c("Outpatient", "Inpatient"), col = c("cornflowerblue", "tan2"), cex = 1, pch = 16, bty = "n")
for (pid in patient_ids) {
  plotColor<-ifelse(metaData$ptInOut[which(metaData$ID == pid)] == "Inpatient", "tan2", "cornflowerblue")
  patient_samples <- brackenM[grep(pid, rownames(brackenM)), grep(pid, colnames(brackenM))]
  pre_sample <- patient_samples[grep("PRE", rownames(patient_samples)), , drop = F]
  distance<-pre_sample[, -grep("PRE", colnames(pre_sample)), drop = F]
  colnames(distance)<-sapply(str_split(colnames(distance), "_", n = 2), `[`, 2)
  colOrder<-natural_sort(colnames(distance))
  distance<-distance[, colOrder, drop = F]
  
  x_positions <- match(colnames(distance), natural_sort(bins))
  y_values <- as.numeric(distance[1, ])
  
  points(x_positions, y_values, pch = 19, col = plotColor)
  
  if(length(y_values) > 1) {
    lines(x_positions, y_values, col = plotColor)
  }

}

for (pid in patient_ids) {
  plotColor<-ifelse(metaData$ptInOut[which(metaData$ID == pid)] == "Inpatient", "tan2", "cornflowerblue")
  patient_samples <- brackenM[grep(pid, rownames(brackenM)), grep(pid, colnames(brackenM))]
  pre_sample <- patient_samples[grep("PRE", rownames(patient_samples)), , drop = F]
  distance<-pre_sample[, -grep("PRE", colnames(pre_sample)), drop = F]
  colnames(distance)<-sapply(str_split(colnames(distance), "_", n = 2), `[`, 2)
  colOrder<-natural_sort(colnames(distance))
  distance<-distance[, colOrder, drop = F]
  
  plot(1, type = "n", xlim = c(1, length(bins)), ylim = range(brackenM),
       xlab = "Date Bin", ylab = "Bray-Curtis Distance", axes = FALSE,
       main = pid)
  axis(1, at = 1:length(bins), labels = natural_sort(bins))
  axis(2)  # Add y-axis
  box()
  
  legend("topright", c("Outpatient", "Inpatient"), col = c("cornflowerblue", "tan2"), cex = 1, pch = 16, bty = "n")
  
  x_positions <- match(colnames(distance), natural_sort(bins))
  y_values <- as.numeric(distance[1, ])
  
  points(x_positions, y_values, pch = 19, col = plotColor)
  
  if(length(y_values) > 1) {
    lines(x_positions, y_values, col = plotColor)
  }
  
}

dev.off()

