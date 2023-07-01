#microbiota shifts between days and baseline 
rm(list = ls())

metaData<-read.csv("metaWithBins.csv", header = T, row.names = 1)

brackenM<-read.csv("distanceMatrix/brackenDistanceMatrix.csv", row.names = 1, check.names = F)

patients<-unique(metaData$ID)

pres<-vector()
posts<-vector()
for (i in 1:length(patients)) {
  patient<-metaData[metaData$ID == patients[i], , drop = F]
  pres[i]<-rownames(patient)[patient$bins == "PRE"]
  posts[i]<-rownames(patient)[which.min(abs(patient$Timepoint - 28))]
}

distance<-vector()
inout<-vector()
for (i in 1:length(patients)) {
  distance[i]<-brackenM[pres[i], posts[i]]
  inout[i]<-metaData[pres[i], 7]
}

names(distance)<-inout

pdf("Plots/PreToD28BrayCurtisDistanceBoxPlots.pdf", width=7, height=7)
par(mar=c(5,6,4,1)+.1)
boxplot(list(Inpatient = distance[names(distance) == "Inpatient"], Outpatient = distance[names(distance) == "Outpatient"]), 
        outline = F,
        cex.lab = 1.5, cex.axis = 1.4, cex.main = 1.5, 
        col=c("olivedrab3", "cornflowerblue"),
        ylab = "Distance between PRE to D28",
        xlab = "Species",
        main = paste0("(Wilcox)P-value=", 
                      signif(wilcox.test(distance[names(distance) == "Inpatient"], distance[names(distance) == "Outpatient"])$p.value, 4)))
stripchart(list(Inpatient = distance[names(distance) == "Inpatient"], Outpatient = distance[names(distance) == "Outpatient"]), 
           method="jitter", 
           vertical=T, pch=19, add=T)

dev.off()