rm(list = ls())

metaData<-read.csv("metaWithBins.csv", header = T, row.names = 1)

#counts tables
brackenT<-read.csv("CountsTables/brackenFiltered.csv", header = T, row.names = 1, check.names = F)
amrT<-read.csv("CountsTables/amrFiltered.csv", header = T, row.names = 1, check.names = F)
rgiT<-read.csv("CountsTables/rgiFiltered.csv", header = T, row.names = 1, check.names = F)
vsearchT<-read.csv("CountsTables/vsearchFiltered.csv", header = T, row.names = 1, check.names = F)

metaBRACKEN<-metaData[colnames(brackenT), , drop = F]
metaAMR<-metaData[colnames(amrT), , drop = F]
metaRGI<-metaData[colnames(rgiT), , drop = F]
metaVSEARCH<-metaData[colnames(vsearchT), , drop = F]

binOrder<-c("PRE", "D0", "D7", "D14", "D21", "D28", "D35", "D60", "D100", "D180")

pdf("Plots/TotalCountsByTimePointsBins.pdf", width=12, height=12)
par(mfrow=c(2, 2))
par(mar=c(5, 6, 4, 1)+.1)

binCounts<-aggregate(colSums(brackenT), list(metaBRACKEN[colnames(brackenT), 6]), FUN=sum)
binCounts$Group.1 <- factor(binCounts$Group.1, levels = binOrder)
binCounts<-binCounts[order(binCounts$Group.1), ]
barplot(binCounts$x, names.arg=c("PRE", "D0", "D7", "D14", "D21", "D28", "D35", "D60", "D100", "D180"), 
        main = "Taxa Counts by Time (Bracken-Species)", cex.axis = 1.5, cex.lab = 1.5, cex.names = 0.8, cex.main = 1.8,
        xlab = "Time Points (Days)", ylab = "Counts(Log10)", col = "tan2")

binCounts<-aggregate(colSums(amrT), list(metaAMR[colnames(amrT), 6]), FUN=sum)
binCounts$Group.1 <- factor(binCounts$Group.1, levels = binOrder)
binCounts<-binCounts[order(binCounts$Group.1), ]
barplot(binCounts$x, names.arg=c("PRE", "D0", "D7", "D14", "D21", "D28", "D35", "D60", "D100", "D180"), 
        main = "Gene Counts by Time (AMR)", cex.axis = 1.5, cex.lab = 1.5, cex.names = 0.8, cex.main = 1.8,
        xlab = "Time Points (Days)", ylab = "Counts(Log10)", col = "coral3")

binCounts<-aggregate(colSums(rgiT), list(metaRGI[colnames(rgiT), 6]), FUN=sum)
binCounts$Group.1 <- factor(binCounts$Group.1, levels = binOrder)
binCounts<-binCounts[order(binCounts$Group.1), ]
barplot(binCounts$x, names.arg=c("PRE", "D0", "D7", "D14", "D21", "D28", "D35", "D60", "D100", "D180"), 
        main = "Gene Counts by Time (RGI)", cex.axis = 1.5, cex.lab = 1.5, cex.names = 0.8, cex.main = 1.8,
        xlab = "Time Points (Days)", ylab = "Counts(Log10)", col = "cornflowerblue")

binCounts<-aggregate(colSums(vsearchT), list(metaVSEARCH[colnames(vsearchT), 6]), FUN=sum)
binCounts$Group.1 <- factor(binCounts$Group.1, levels = binOrder)
binCounts<-binCounts[order(binCounts$Group.1), ]
barplot(binCounts$x, names.arg=c("PRE", "D0", "D7", "D14", "D21", "D28", "D35", "D60", "D100", "D180"), 
        main = "Gene Counts by Time (vsearch)", cex.axis = 1.5, cex.lab = 1.5, cex.names = 0.8, cex.main = 1.8,
        xlab = "Time Points (Days)", ylab = "Counts(Log10)", col = "olivedrab4")

dev.off()
