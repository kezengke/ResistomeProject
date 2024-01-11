#bray-curtis distance from pre to other bin for each patient (normalized, no log for plots)
rm(list = ls())
library("vegan")
library(dplyr)
library(ggplot2)
library(stringr)

metaData<-read.csv("MetaWithInOutDonerTreatmentType.csv", header = T, row.names = 1)
metaData$days<-metaData$Timepoint
metaData$days[which(metaData$bins == "PRE")]<-"PRE"

#counts tables
brackenT<-read.csv("CountsTables/brackenProcessed.csv", header = T, row.names = 1, check.names = F)
brackenT<-brackenT[, intersect(rownames(metaData), colnames(brackenT)), drop = F]

Norm <- function(counts) {
  #Normalization
  n<-colSums(counts)
  sumx<-sum(counts)
  for (i in 1:ncol(counts)) {
    counts[,i]<-counts[,i]/n[i]
  }
  counts<-(counts*(sumx/ncol(counts)))
  
  counts<-data.frame(counts, check.names = F)
  return(counts)
}
brackenT<-Norm(brackenT)

distanceMatrix<-function(countsT, dissimilarityMetric){
  distanceT = as.matrix(vegdist(t(countsT), method = dissimilarityMetric))
  return(distanceT)
}

natural_sort <- function(x) {
  numeric_part <- as.numeric(x)
  ordered_indices <- order(numeric_part)
  x[ordered_indices]
}

brackenM<-distanceMatrix(brackenT, "bray")
names<-paste0(metaData[colnames(brackenM), 5], "_",metaData[colnames(brackenM), 11])
colnames(brackenM)<-names
rownames(brackenM) <- make.names(names, unique = FALSE)

patient_ids<-unique(metaData$ID)

pdf("Plots/(daysNoLog)BrayCurtisDistanceFromBaselineForEachPatient.pdf", width=20, height=15)
par(mfrow=c(3, 2))
par(mar=c(5, 6, 4, 1)+.1)

nonPreSamples<-metaData[which(metaData$bins != "PRE"), , drop = F]

plot(1, type = "n", xlim = range(nonPreSamples$Timepoint), ylim = range(brackenM),
     xlab = "Days from PRE", ylab = "Bray-Curtis Distance", axes = FALSE)
axis(1)
axis(2)  # Add y-axis
box()
legend("bottomright", c("Outpatient", "Inpatient"), col = c("cornflowerblue", "tan2"), cex = 1, pch = 16, bty = "n")
for (pid in patient_ids) {
  plotColor<-ifelse(metaData$ptInOut[which(metaData$ID == pid)] == "Inpatient", "tan2", "cornflowerblue")
  patient_samples <- brackenM[grep(pid, rownames(brackenM)), grep(pid, colnames(brackenM))]
  pre_sample <- patient_samples[grep("PRE", rownames(patient_samples)), , drop = F]
  distance<-pre_sample[, -grep("PRE", colnames(pre_sample)), drop = F]
  colnames(distance)<-sapply(str_split(colnames(distance), "_", n = 2), `[`, 2)
  colOrder<-natural_sort(colnames(distance))
  distance<-distance[, colOrder, drop = F]
  colnames(distance)<-as.numeric(colnames(distance))
  
  points(as.numeric(colnames(distance)), as.numeric(distance[1, ]), pch = 19, col = plotColor)
  
  if(length(as.numeric(distance[1, ])) > 1) {
    lines(as.numeric(colnames(distance)), as.numeric(distance[1, ]), col = plotColor)
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
  colnames(distance)<-as.numeric(colnames(distance))
  
  plot(1, type = "n", xlim = range(nonPreSamples$Timepoint), ylim = range(brackenM),
       xlab = "Days from PRE", ylab = "Bray-Curtis Distance", 
       main = pid, axes = FALSE)
  axis(1)
  axis(2)  # Add y-axis
  box()
  
  legend("bottomright", c("Outpatient", "Inpatient"), col = c("cornflowerblue", "tan2"), cex = 1, pch = 16, bty = "n")
  
  points(as.numeric(colnames(distance)), as.numeric(distance[1, ]), pch = 19, col = plotColor)
  
  if(length(as.numeric(distance[1, ])) > 1) {
    lines(as.numeric(colnames(distance)), as.numeric(distance[1, ]), col = plotColor)
  }
  
}

dev.off()

