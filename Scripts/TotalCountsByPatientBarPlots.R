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

pdf("Plots/TotalCountsByPatients.pdf", width=12, height=12)
par(mfrow=c(2, 2))
par(mar=c(5, 6, 4, 1)+.1)

binCounts<-aggregate(colSums(brackenT), list(metaBRACKEN[colnames(brackenT), 5]), FUN=sum)
barplot(binCounts$x, names.arg=binCounts$Group.1, las = 3,
        main = "Taxa Counts by Patient (Bracken-Species)", col = "tan2",
        ylab = "Taxa Counts(Log10)",
        cex.axis = 1, cex.lab = 1.5, cex.names = 1, cex.main = 1.8)

binCounts<-aggregate(colSums(amrT), list(metaAMR[colnames(amrT), 5]), FUN=sum)
barplot(binCounts$x, names.arg=binCounts$Group.1, las = 3, 
        main = "Gene Counts by Patient (AMR)", col = "coral3", 
        ylab = "Gene Counts(Log10)",
        cex.axis = 1, cex.lab = 1.5, cex.names = 1, cex.main = 1.8)

binCounts<-aggregate(colSums(rgiT), list(metaRGI[colnames(rgiT), 5]), FUN=sum)
barplot(binCounts$x, names.arg=binCounts$Group.1, las = 3, 
        main = "Gene Counts by Patient (RGI)", col = "cornflowerblue",
        ylab = "Gene Counts(Log10)", 
        cex.axis = 1, cex.lab = 1.5, cex.names = 1, cex.main = 1.8)

binCounts<-aggregate(colSums(vsearchT), list(metaVSEARCH[colnames(vsearchT), 5]), FUN=sum)
barplot(binCounts$x, names.arg=binCounts$Group.1, las = 3,
        main = "Gene Counts by Patient (vsearch)", col = "olivedrab4",
        ylab = "Gene Counts(Log10)", 
        cex.axis = 1, cex.lab = 1.5, cex.names = 1, cex.main = 1.8,
)

dev.off()
