#pcoa for in vs out patient for pre and D28 samples
rm(list = ls())
library("vegan")

metaData<-read.csv("metaWithBins.csv", header = T, row.names = 1)

#counts tables
brackenT<-read.csv("CountsTables/brackenProcessed.csv", header = T, row.names = 1, check.names = F)

patients<-unique(metaData$ID)
pres<-vector()
posts<-vector()
for (i in 1:length(patients)) {
  patient<-metaData[metaData$ID == patients[i], , drop = F]
  pres[i]<-rownames(patient)[patient$bins == "PRE"]
  posts[i]<-rownames(patient)[which.min(abs(patient$Timepoint - 28))]
}

brackenPre<-brackenT[, pres, drop = F]
brackenPost<-brackenT[, posts, drop = F]
metaPre<-metaData[pres, , drop = F]
metaPost<-metaData[posts, , drop = F]

pdf("Plots/MDSForPreAndD28InAndOutPatients.pdf", width=12, height=6)
par(mfrow=c(1,2))
par(mar=c(5,6,4,1)+.1)
circleCol<-c("cornflowerblue", "darkorange")
cols<-circleCol[factor(metaPre$ptInOut, levels = c("Outpatient", "Inpatient"))]
MDS<-capscale(t(brackenPre)~1,distance = "bray")
percentVariance<-MDS$CA$eig/sum(eigenvals(MDS))*100
pval<-adonis2(t(brackenPre)~metaPre$ptInOut, method="bray")$`Pr(>F)`[1]
statusPlot<-ordiplot(MDS,choices = c(1,2),type="none",cex.lab = 1,
                     xlab = paste("MDS1  ", format(percentVariance[1], digits = 4), "%", sep = ""),
                     ylab = paste("MDS2  ", format(percentVariance[2], digits = 4), "%", sep = ""),
                     main = paste("Species-PRE n =", ncol(brackenPre), "\nP-value:",pval))
points(statusPlot,"sites", pch = 19, cex = 2.5, col = adjustcolor(cols, alpha.f = 0.5))
ordiellipse(statusPlot, metaPre$ptInOut, kind="se", conf=0.95, lwd=4, draw = "lines", col=circleCol[1],show.groups="Outpatient",label=T,font=2,cex=1) 
ordiellipse(statusPlot, metaPre$ptInOut, kind="se", conf=0.95, lwd=4, draw = "lines", col=circleCol[2],show.groups="Inpatient",label=T,font=2,cex=1) 
legend("topright", c("Outpatient", "Inpatient"), col = circleCol[1:2], cex = 1, pch = 16, bty = "n")

circleCol<-c("cornflowerblue", "darkorange")
cols<-circleCol[factor(metaPost$ptInOut, levels = c("Outpatient", "Inpatient"))]
MDS<-capscale(t(brackenPost)~1,distance = "bray")
percentVariance<-MDS$CA$eig/sum(eigenvals(MDS))*100
pval<-adonis2(t(brackenPost)~metaPost$ptInOut, method="bray")$`Pr(>F)`[1]
statusPlot<-ordiplot(MDS,choices = c(1,2),type="none",cex.lab = 1,
                     xlab = paste("MDS1  ", format(percentVariance[1], digits = 4), "%", sep = ""),
                     ylab = paste("MDS2  ", format(percentVariance[2], digits = 4), "%", sep = ""),
                     main = paste("Species-D28 n =", ncol(brackenPost), "\nP-value:",pval))
points(statusPlot,"sites", pch = 19, cex = 2.5, col = adjustcolor(cols, alpha.f = 0.5))
ordiellipse(statusPlot, metaPost$ptInOut, kind="se", conf=0.95, lwd=4, draw = "lines", col=circleCol[1],show.groups="Outpatient",label=T,font=2,cex=1) 
ordiellipse(statusPlot, metaPost$ptInOut, kind="se", conf=0.95, lwd=4, draw = "lines", col=circleCol[2],show.groups="Inpatient",label=T,font=2,cex=1) 
legend("topright", c("Outpatient", "Inpatient"), col = circleCol[1:2], cex = 1, pch = 16, bty = "n")

dev.off()