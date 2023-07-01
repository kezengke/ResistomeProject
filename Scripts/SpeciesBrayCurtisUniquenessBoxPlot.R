#uniqueness boxplots
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

timePoint<-c(rep("PRE", length(pres)), rep("D28", length(posts)))
timePoint<-data.frame(timePoint)
rownames(timePoint)<-c(pres, posts)
timePoint$InOut<-metaData[rownames(timePoint), 7]
newMeta<-paste(timePoint$timePoint, timePoint$InOut, sep = "")
newMeta<-data.frame(newMeta)
rownames(newMeta)<-rownames(timePoint)
newMeta<-newMeta[order(newMeta$newMeta), , drop = F]

smallest_non_zero<-function(x) {
  x <- x[x != 0]  #exclude zeros
  if(length(x) > 0) min(x) else NA# return smallest non-zero value if it exists, otherwise return NA
}

closestDistance<-apply(brackenM, 1, smallest_non_zero)
closestDistance<-data.frame(closestDistance)
rownames(closestDistance)<-rownames(brackenM)
closestDistance<-closestDistance[rownames(newMeta), , drop = F]

pdf("Plots/UniquenessBrayCurtisDistanceBoxPlots.pdf", width=7, height=7)
par(mar=c(5,6,4,1)+.1)
boxplot(unlist(closestDistance)~newMeta$newMeta, 
        outline = F, names = c("D28In", "D28Out", "PREIn", "PREOut"),
        cex.lab = 1.5, cex.axis = 1, cex.main = 1.5, 
        col=c("cornflowerblue", "darkorange", "coral3", "olivedrab3"),
        ylab = "Bray-Curtis Uniqueness",
        xlab = "Patient",
        main = paste0("Species(ANOVA)P-value=", 
                      signif(summary(aov(unlist(closestDistance)~newMeta$newMeta))[[1]]["newMeta$newMeta", "Pr(>F)"], 3)))
stripchart(unlist(closestDistance)~newMeta$newMeta, 
           method="jitter", 
           vertical=T, pch=19, add=T)

dev.off()

