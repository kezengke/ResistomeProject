#barplots for each patient and their taxa/gene type richness
rm(list = ls())

metaData<-read.csv("metaWithBins.csv", header = T, row.names = 1)

pdf("Plots/TypeRichnessByPatients.pdf", width=12, height=12)
par(mfrow=c(2, 2))
par(mar=c(5, 6, 4, 1)+.1)
#BRACKEN
brackenT<-read.csv("CountsTables/brackenFiltered.csv", header = T, row.names = 1, check.names = F)
metaBRACKEN<-metaData[colnames(brackenT), , drop = F]

binCounts<-aggregate(colSums(brackenT>0), list(metaBRACKEN[colnames(brackenT), 5]), FUN=sum)
barplot(binCounts$x, names.arg=binCounts$Group.1, las = 3,
        main = "Taxa Richness by Patient (Bracken-Species)", col = adjustcolor("tan2", alpha.f = 0.7),
        ylab = "Type of Taxa",
        cex.axis = 1, cex.lab = 1.5, cex.names = 1, cex.main = 1.7)

#AMR
amrT<-read.csv("CountsTables/amrFiltered.csv", header = T, row.names = 1, check.names = F)
metaAMR<-metaData[colnames(amrT), , drop = F]

binCounts<-aggregate(colSums(amrT>0), list(metaAMR[colnames(amrT), 5]), FUN=sum)
barplot(binCounts$x, names.arg=binCounts$Group.1, las = 3,
        main = "Gene Richness by Patient (AMR)", col = adjustcolor("coral3", alpha.f = 0.7),
        ylab = "Type of Genes",
        cex.axis = 1, cex.lab = 1.5, cex.names = 1, cex.main = 1.8)

#RGI
rgiT<-read.csv("CountsTables/rgiFiltered.csv", header = T, row.names = 1, check.names = F)
metaRGI<-metaData[colnames(rgiT), , drop = F]

binCounts<-aggregate(colSums(rgiT>0), list(metaRGI[colnames(rgiT), 5]), FUN=sum)
barplot(binCounts$x, names.arg=binCounts$Group.1, las = 3,
        main = "Gene Richness by Patient (RGI)", col = adjustcolor("cornflowerblue", alpha.f = 0.7),
        ylab = "Type of Genes",
        cex.axis = 1, cex.lab = 1.5, cex.names = 1, cex.main = 1.8)

#VSEARCH
vsearchT<-read.csv("CountsTables/vsearchFiltered.csv", header = T, row.names = 1, check.names = F)
metaVSEARCH<-metaData[colnames(vsearchT), , drop = F]

binCounts<-aggregate(colSums(vsearchT>0), list(metaVSEARCH[colnames(vsearchT), 5]), FUN=sum)
barplot(binCounts$x, names.arg=binCounts$Group.1, las = 3,
        main = "Gene Richness by Patient (VSEARCH)", col = adjustcolor("olivedrab4", alpha.f = 0.7),
        ylab = "Type of Genes",
        cex.axis = 1, cex.lab = 1.5, cex.names = 1, cex.main = 1.8)

dev.off()
