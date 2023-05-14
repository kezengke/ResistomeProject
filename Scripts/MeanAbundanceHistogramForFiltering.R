pdf("Plots/MeanAbundance(Log).pdf", width=21, height=14)
par(mfrow=c(2,3))
#Bracken
rm(list = ls())
myT<-read.csv("CountsTables/brackenNormalized.csv", row.name = 1, header = T, check.names = F)
hist.data<-hist(rowMeans(myT), breaks = 50, plot=F)
hist.data$counts<-log10(hist.data$counts+1)
plot(hist.data, cex.lab=1.5,
     main = "Bracken(species)", 
     xlab = "taxa mean abundance", ylab = "log10(Frequency)")
#AMR
rm(list = ls())
myT<-read.csv("CountsTables/amrNormalized.csv", row.name = 1, header = T, check.names = F)
hist.data<-hist(rowMeans(myT), breaks = 50, plot=F)
hist.data$counts<-log10(hist.data$counts+1)
plot(hist.data, cex.lab=1.5,
     main = "AMR", 
     xlab = "gene mean abundance", ylab = "log10(Frequency)")
#RGI
rm(list = ls())
myT<-read.csv("CountsTables/rgiNormalized.csv", row.name = 1, header = T, check.names = F)
hist.data<-hist(rowMeans(myT), breaks = 50, plot=F)
hist.data$counts<-log10(hist.data$counts+1)
plot(hist.data, cex.lab=1.5,
     main = "RGI", 
     xlab = "gene mean abundance", ylab = "log10(Frequency)")

#VSEARCH
rm(list = ls())
myT<-read.csv("CountsTables/vsearchNormalized.csv", header = T, row.names = 1, check.names = F)
hist.data<-hist(rowMeans(myT), breaks = 50, plot=F)
hist.data$counts<-log10(hist.data$counts+1)
plot(hist.data, cex.lab=1.5,
     main = "vsearch", 
     xlab = "gene mean abundance", ylab = "log10(Frequency)")

dev.off()

