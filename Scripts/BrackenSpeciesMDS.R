rm(list = ls())
library(stringr)
brackenT<-read.delim("CountsTables/duke_bracken.csv", sep = ",", header = T, row.names = 1)
brackenT<-brackenT[, -c(1,2)] #get rid of taxonomy ID and level
brackenT<-brackenT[, grepl("num", colnames(brackenT))] #get rid of fractions

#Normalization
n<-colSums(brackenT)
sumx<-sum(brackenT)
brackenT<-log10((brackenT/n)*(sumx/ncol(brackenT))+1)
brackenT<-data.frame(brackenT)

library("vegan")

pdf("Plots/BrackenSpeciesMDS.pdf",onefile = T)
par(mfrow=c(1,1))
MDS<-capscale(t(brackenT)~1,distance = "bray")
percentVariance<-MDS$CA$eig/sum(eigenvals(MDS))*100
statusPlot<-ordiplot(MDS,choices=c(1,2),type="none",cex.lab=1,
                     xlab=paste("MDS1  ", format(percentVariance[1], digits = 4), "%", sep = ""),
                     ylab=paste("MDS2  ", format(percentVariance[2], digits = 4), "%", sep = ""),
                     main="Bracken-Species ")
points(statusPlot,"sites", pch=19, cex=2.5, col=adjustcolor("cornflowerblue", alpha.f = 0.5))

dev.off()