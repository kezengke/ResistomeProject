rm(list = ls())

metaData<-read.csv("metaWithBins.csv", header = T, row.names = 1)

#counts tables
brackenT<-read.csv("CountsTables/brackenNormalized.csv", header = T, row.names = 1, check.names = F)
amrT<-read.csv("CountsTables/amrNormalized.csv", header = T, row.names = 1, check.names = F)
rgiT<-read.csv("CountsTables/rgiNormalized.csv", header = T, row.names = 1, check.names = F)
vsearchT<-read.csv("CountsTables/vsearchNormalized.csv", header = T, row.names = 1, check.names = F)

metaBRACKEN<-metaData[colnames(brackenT), , drop = F]
metaAMR<-metaData[colnames(amrT), , drop = F]
metaRGI<-metaData[colnames(rgiT), , drop = F]
metaVSEARCH<-metaData[colnames(vsearchT), , drop = F]

pdf("Plots/TotalCountsByTimePointsBins.pdf", width=12, height=12)
par(mfrow=c(2, 2))
par(mar=c(5, 6, 4, 1)+.1)

binCounts<-aggregate(colSums(brackenT), list(metaBRACKEN[colnames(brackenT), 5]), FUN=sum)
binCounts<-binCounts[c(7, 6, 2, 3, 4, 5, 1),]
barplot(binCounts$x, names.arg=c("PRE", "D7", "D14", "D21", "D35", "D60", "D100"), 
        main = "Taxa Counts by Time (Bracken-Species)", cex.axis = 1.5, cex.lab = 1.5, cex.names = 0.8, cex.main = 1.8,
        xlab = "Time Points (Days)", ylab = "Counts(Log10)", col = "tan2")

binCounts<-aggregate(colSums(amrT), list(metaAMR[colnames(amrT), 5]), FUN=sum)
binCounts<-binCounts[c(7, 6, 2, 3, 4, 5, 1),]
barplot(binCounts$x, names.arg=c("PRE", "D7", "D14", "D21", "D35", "D60", "D100"), 
        main = "Gene Counts by Time (AMR)", cex.axis = 1.5, cex.lab = 1.5, cex.names = 0.8, cex.main = 1.8,
        xlab = "Time Points (Days)", ylab = "Counts(Log10)", col = "coral3")

binCounts<-aggregate(colSums(rgiT), list(metaRGI[colnames(rgiT), 5]), FUN=sum)
binCounts<-binCounts[c(7, 6, 2, 3, 4, 5, 1),]
barplot(binCounts$x, names.arg=c("PRE", "D7", "D14", "D21", "D35", "D60", "D100"), 
        main = "Gene Counts by Time (RGI)", cex.axis = 1.5, cex.lab = 1.5, cex.names = 0.8, cex.main = 1.8,
        xlab = "Time Points (Days)", ylab = "Counts(Log10)", col = "cornflowerblue")

binCounts<-aggregate(colSums(vsearchT), list(metaVSEARCH[colnames(vsearchT), 5]), FUN=sum)
binCounts<-binCounts[c(7, 6, 2, 3, 4, 5, 1),]
barplot(binCounts$x, names.arg=c("PRE", "D7", "D14", "D21", "D35", "D60", "D100"), 
        main = "Gene Counts by Time (vsearch)", cex.axis = 1.5, cex.lab = 1.5, cex.names = 0.8, cex.main = 1.8,
        xlab = "Time Points (Days)", ylab = "Counts(Log10)", col = "olivedrab4")

dev.off()
