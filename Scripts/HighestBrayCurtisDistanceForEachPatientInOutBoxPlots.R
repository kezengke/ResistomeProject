#highest bray-curtis distance from pre box plots for in and out patients
rm(list = ls())
library("vegan")
library(ggplot2)
library(stringr)
library(gridExtra)

metaData<-read.csv("MetaWithInOutDonerTreatmentType.csv", header = T, row.names = 1)
metaData$days<-metaData$Timepoint
metaData$days[which(metaData$bins == "PRE")]<-"PRE"
#remove duplicated pre samples
metaData<-metaData[metaData$days > 0, , drop = F]

# #keep only Allo samples
# metaData<-metaData[which(metaData$donor == "Allo"), , drop = F]

#calculate bray-curtis distance
distanceMatrix<-function(countsT, dissimilarityMetric){
  distanceT = as.matrix(vegdist(t(countsT), method = dissimilarityMetric))
  return(distanceT)
}

#find the max distance among each patient
distProcessor <- function(distMat, meta) {
  names<-paste0(meta[colnames(distMat), 5], "_",meta[colnames(distMat), 12])
  colnames(distMat)<-names
  rownames(distMat) <- make.names(names, unique = FALSE)
  
  patient_ids<-unique(meta$ID)
  max_distance<-vector()
  for (pid in patient_ids) {
    #find all the samples for the same patient
    patient_samples <- distMat[grep(pid, rownames(distMat)), grep(pid, colnames(distMat))]
    #locate the distance from pre
    pre_sample <- patient_samples[grep("PRE", rownames(patient_samples)), , drop = F]
    #exclude pre to pre
    distance<-pre_sample[, -grep("PRE", colnames(pre_sample)), drop = F]
    
    max_distance[pid]<-max(distance)
  }
  return(max_distance)
}

plotFun <- function(distMat, meta, type, tableType) {
  plotMeta<-unique(data.frame(meta$ID, meta[[type]]))
  rownames(plotMeta)<-plotMeta[,1]
  colnames(plotMeta)<-c("ID", "type")
  pval1<-t.test(distMat~plotMeta$type)$p.value
  pval2<-wilcox.test(distMat~plotMeta$type)$p.value
  plotMeta<-plotMeta[, 2, drop = F]
  
  p<-ggplot(plotMeta, aes(x = type, y = distMat, fill = type)) +
    geom_boxplot() +  
    geom_jitter(width = 0.2, color = "black", size = 1.5) +  
    scale_fill_manual(values = c("olivedrab3", "cornflowerblue")) + 
    theme_classic() +
    labs(
      x = "Patient type",
      y = "Max distance per patient",
      title = paste0(tableType, "\nMax bray-curtis distance for each patient \nT-test pval=", signif(pval1, 4),
                     "\nWilcox pval=", signif(pval2, 4))
    ) +
    theme(
      legend.position="none",
      plot.title = element_text(size = 13, hjust = 0.5), # Adjust title size and position
      axis.title.x = element_text(size = 14), # Adjust x axis label size
      axis.title.y = element_text(size = 14), # Adjust y axis label size
      axis.text.x = element_text(size = 14), # Adjust x axis tick label size
      axis.text.y = element_text(size = 14) # Adjust y axis tick label size
    )
  
  return(p)
}

pdf("Plots/(Filtered)MaxBrayCurtisDistanceForEachPatientBoxPlots.pdf", width=18, height=6)
par(mar=c(5,6,4,1)+.1)

#species
speciesT<-read.csv("CountsTables/brackenFiltered.csv", header = T, row.names = 1, check.names = F)
speciesT<-speciesT[, intersect(rownames(metaData), colnames(speciesT)), drop = F]

myM<-distanceMatrix(speciesT, "bray")
myDist<-distProcessor(myM, metaData)
p1<-plotFun(myDist, metaData, "ptInOut", "Species")

#genus
genusT<-read.csv("CountsTables/genusFiltered.csv", header = T, row.names = 1, check.names = F)
genusT<-genusT[, intersect(rownames(metaData), colnames(genusT)), drop = F]

myM<-distanceMatrix(genusT, "bray")
myDist<-distProcessor(myM, metaData)
p2<-plotFun(myDist, metaData, "ptInOut", "Genus")

#phylum
phylumT<-read.csv("CountsTables/phylumFiltered.csv", header = T, row.names = 1, check.names = F)
phylumT<-phylumT[, intersect(rownames(metaData), colnames(phylumT)), drop = F]

myM<-distanceMatrix(phylumT, "bray")
myDist<-distProcessor(myM, metaData)
p3<-plotFun(myDist, metaData, "ptInOut", "Phylum")

grid.arrange(p1, p2, p3, nrow = 1)

#amr
amrT<-read.csv("CountsTables/amrFiltered.csv", header = T, row.names = 1, check.names = F)
amrT<-amrT[, intersect(rownames(metaData), colnames(amrT)), drop = F]

myM<-distanceMatrix(amrT, "bray")
myDist<-distProcessor(myM, metaData)
p4<-plotFun(myDist, metaData, "ptInOut", "AMR")

#rgi
rgiT<-read.csv("CountsTables/rgiFiltered.csv", header = T, row.names = 1, check.names = F)
rgiT<-rgiT[, intersect(rownames(metaData), colnames(rgiT)), drop = F]

myM<-distanceMatrix(rgiT, "bray")
myDist<-distProcessor(myM, metaData)
p5<-plotFun(myDist, metaData, "ptInOut", "RGI")

#vsearch
vsearchT<-read.csv("CountsTables/vsearchFiltered.csv", header = T, row.names = 1, check.names = F)
vsearchT<-vsearchT[, intersect(rownames(metaData), colnames(vsearchT)), drop = F]

myM<-distanceMatrix(vsearchT, "bray")
myDist<-distProcessor(myM, metaData)
p6<-plotFun(myDist, metaData, "ptInOut", "Vsearch")

grid.arrange(p4, p5, p6, nrow = 1)

# #pathway
# pathwayT<-read.csv("CountsTables/pathFiltered.csv", header = T, row.names = 1, check.names = F)
# pathwayT<-pathwayT[, intersect(rownames(metaData), colnames(pathwayT)), drop = F]
# 
# myM<-distanceMatrix(pathwayT, "bray")
# myDist<-distProcessor(myM, metaData)
# p7<-plotFun(myDist, metaData, "ptInOut", "Pathway")
# 
# grid.arrange(p7, nrow = 1, ncol = 3)

dev.off()
