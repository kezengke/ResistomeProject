#pcoa for in vs out patient for pre and D28 samples
rm(list = ls())
library("vegan")

metaData<-read.csv("metaWithBins.csv", header = T, row.names = 1)

#counts tables
rgiT<-read.csv("CountsTables/rgiProcessed.csv", header = T, row.names = 1, check.names = F)

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

myT<-rgiT[, rownames(newMeta), drop = F]

pdf("Plots/MDSForPreAndD28InAndOutPatientsRGI(FourColors).pdf", width=6, height=6)
par(mfrow=c(1,1))
par(mar=c(5,6,4,1)+.1)
circleCol<-c("cornflowerblue", "darkorange", "coral3", "olivedrab3")
cols<-circleCol[factor(newMeta$newMeta, levels = c("D28Inpatient", "D28Outpatient", "PREInpatient", "PREOutpatient"))]
MDS<-capscale(t(myT)~1,distance = "bray")
percentVariance<-MDS$CA$eig/sum(eigenvals(MDS))*100
pval<-adonis2(t(myT)~newMeta$newMeta, method="bray")$`Pr(>F)`[1]
statusPlot<-ordiplot(MDS,choices = c(1,2),type="none",cex.lab = 1,
                     xlab = paste("MDS1  ", format(percentVariance[1], digits = 4), "%", sep = ""),
                     ylab = paste("MDS2  ", format(percentVariance[2], digits = 4), "%", sep = ""),
                     main = paste("RGI n =", ncol(myT), "\nP-value:",pval))
points(statusPlot,"sites", pch = 19, cex = 2.5, col = adjustcolor(cols, alpha.f = 0.5))
ordiellipse(statusPlot, newMeta$newMeta, kind="se", conf=0.95, lwd=4, draw = "lines", col=circleCol[1],show.groups="D28Inpatient",label=T,font=2,cex=1) 
ordiellipse(statusPlot, newMeta$newMeta, kind="se", conf=0.95, lwd=4, draw = "lines", col=circleCol[2],show.groups="D28Outpatient",label=T,font=2,cex=1)
ordiellipse(statusPlot, newMeta$newMeta, kind="se", conf=0.95, lwd=4, draw = "lines", col=circleCol[3],show.groups="PREInpatient",label=T,font=2,cex=1) 
ordiellipse(statusPlot, newMeta$newMeta, kind="se", conf=0.95, lwd=4, draw = "lines", col=circleCol[4],show.groups="PREOutpatient",label=T,font=2,cex=1) 
legend("topright", c("D28Inpatient", "D28Outpatient", "PREInpatient", "PREOutpatient"), col = circleCol[1:4], cex = 1, pch = 16, bty = "n")

dev.off()