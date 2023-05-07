rm(list = ls())
library("nlme")

metaData<-read.csv("metaWithBins.csv", header = T, row.names = 1)

#counts tables
brackenT<-read.csv("CountsTables/brackenNormalized.csv", header = T, row.names = 1, check.names = F)
amrT<-read.csv("CountsTables/amrNormalized.csv", header = T, row.names = 1, check.names = F)
rgiT<-read.csv("CountsTables/rgiNormalized.csv", header = T, row.names = 1, check.names = F)
vsearchT<-read.csv("CountsTables/vsearchNormalized.csv", header = T, row.names = 1, check.names = F)

#filtering
brackenT<-brackenT[apply(brackenT == 0, 1, sum) <= (ncol(brackenT)*0.8), ]
amrT<-amrT[apply(amrT == 0, 1, sum) <= (ncol(amrT)*0.8), ]
rgiT<-rgiT[apply(rgiT == 0, 1, sum) <= (ncol(rgiT)*0.8), ]
vsearchT<-vsearchT[apply(vsearchT == 0, 1, sum) <= (ncol(vsearchT)*0.8), ]

metaBRACKEN<-metaData[colnames(brackenT), , drop = F]
metaAMR<-metaData[colnames(amrT), , drop = F]
metaRGI<-metaData[colnames(rgiT), , drop = F]
metaVSEARCH<-metaData[colnames(vsearchT), , drop = F]

#lm 2nd order
pdf("Plots/SecondOrderLm(Bracken).pdf", width=12, height=18)
par(mfrow=c(3, 2))
par(mar=c(5, 6, 4, 1)+.1)

for (i in 1:nrow(brackenT)) {
  Model<-lm(unlist(brackenT[i,]) ~ poly(metaBRACKEN$Timepoint, 2) , x=T)
  plot(metaBRACKEN$Timepoint, unlist(brackenT[i,]), 
       pch = 19, col = "tan2", 
       main = paste("Bracken(Species)\n", rownames(brackenT)[i]),
       xlab = "TimePoints", ylab = "Counts(Log10)")
  lines(sort(metaBRACKEN$Timepoint), fitted(Model)[order(metaBRACKEN$Timepoint)], 
        col = "grey", type = "l")
}
dev.off()

pdf("Plots/SecondOrderLm(AMR).pdf", width=12, height=18)
par(mfrow=c(3, 2))
par(mar=c(5, 6, 4, 1)+.1)
for (i in 1:nrow(amrT)) {
  Model<-lm(unlist(amrT[i,]) ~ poly(metaAMR$Timepoint, 2) , x=T)
  plot(metaAMR$Timepoint, unlist(amrT[i,]), 
       pch = 19, col = "coral3", 
       main = paste("AMR\n", rownames(amrT)[i]),
       xlab = "TimePoints", ylab = "Counts(Log10)")
  lines(sort(metaAMR$Timepoint), fitted(Model)[order(metaAMR$Timepoint)], 
        col = "grey", type = "l")
}
dev.off()

pdf("Plots/SecondOrderLm(RGI).pdf", width=12, height=18)
par(mfrow=c(3, 2))
par(mar=c(5, 6, 4, 1)+.1)
for (i in 1:nrow(rgiT)) {
  Model<-lm(unlist(rgiT[i,]) ~ poly(metaRGI$Timepoint, 2) , x=T)
  plot(metaRGI$Timepoint, unlist(rgiT[i,]), 
       pch = 19, col = "cornflowerblue", 
       main = paste("RGI\n", rownames(rgiT)[i]),
       xlab = "TimePoints", ylab = "Counts(Log10)")
  lines(sort(metaRGI$Timepoint), fitted(Model)[order(metaRGI$Timepoint)], 
        col = "grey", type = "l")
}
dev.off()

pdf("Plots/SecondOrderLm(vsearch).pdf", width=12, height=18)
par(mfrow=c(3, 2))
par(mar=c(5, 6, 4, 1)+.1)
for (i in 1:nrow(vsearchT)) {
  Model<-lm(unlist(vsearchT[i,]) ~ poly(metaVSEARCH$Timepoint, 2) , x=T)
  plot(metaVSEARCH$Timepoint, unlist(vsearchT[i,]), 
       pch = 19, col = "olivedrab4", 
       main = paste("vsearch\n", rownames(vsearchT)[i]),
       xlab = "TimePoints", ylab = "Counts(Log10)")
  lines(sort(metaVSEARCH$Timepoint), fitted(Model)[order(metaVSEARCH$Timepoint)], 
        col = "grey", type = "l")
}
dev.off()