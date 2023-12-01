#in vs out patient wilcoxon test for pre and post patient histograms.
rm(list = ls())

plotFun <- function(resPre, resPost, resType, color) {
        hist(resPre[, 1], breaks=seq(0, 1, 0.05), xlab = "p-value", col = color,
             main = paste0("(", resType, ")", "Pre In vs Out Wilcoxon P-vals"), 
             xlim=c(0,1), cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.4)
        abline(v=0.05, col=gray(.5), lty=2)
        mtext(at=0.05, side=3, text=0.05, col=gray(.5))
        
        hist(resPost[, 1], breaks=seq(0, 1, 0.05), xlab = "p-value", col = color,
             main = paste0("(", resType, ")", "Post In vs Out Wilcoxon P-vals"), 
             xlim=c(0,1), cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.4)
        abline(v=0.05, col=gray(.5), lty=2)
        mtext(at=0.05, side=3, text=0.05, col=gray(.5))
}

#wilcox results
speciesPreResults<-read.csv("Wilcox/speciesPreInOut.csv", row.names = 1)
speciesPostResults<-read.csv("Wilcox/speciesPostInOut.csv", row.names = 1)
genusPreResults<-read.csv("Wilcox/genusPreInOut.csv", row.names = 1)
genusPostResults<-read.csv("Wilcox/genusPostInOut.csv", row.names = 1)
phylumPreResults<-read.csv("Wilcox/phylumPreInOut.csv", row.names = 1)
phylumPostResults<-read.csv("Wilcox/phylumPostInOut.csv", row.names = 1)

amrPreResults<-read.csv("Wilcox/amrPreInOut.csv", row.names = 1)
amrPostResults<-read.csv("Wilcox/amrPostInOut.csv", row.names = 1)
rgiPreResults<-read.csv("Wilcox/rgiPreInOut.csv", row.names = 1)
rgiPostResults<-read.csv("Wilcox/rgiPostInOut.csv", row.names = 1)
vsearchPreResults<-read.csv("Wilcox/vsearchPreInOut.csv", row.names = 1)
vsearchPostResults<-read.csv("Wilcox/vsearchPostInOut.csv", row.names = 1)

pathwayPreResults<-read.csv("Wilcox/pathwayPreInOut.csv", row.names = 1)
pathwayPostResults<-read.csv("Wilcox/pathwayPostInOut.csv", row.names = 1)


pdf("Plots/PrePostInVsOutPatientWilcoxonPvalHistograms.pdf", width=12, height=18)
par(mfrow=c(3, 2))
par(mar=c(5, 6, 4, 1)+.1)

plotFun(speciesPreResults, speciesPostResults, "Species", "tan2")
plotFun(genusPreResults, genusPostResults, "Genus", "orchid3")
plotFun(phylumPreResults, phylumPostResults, "Phylum", "coral3")

plotFun(amrPreResults, amrPostResults, "AMR", "cornflowerblue")
plotFun(rgiPreResults, rgiPostResults, "RGI", "olivedrab4")
plotFun(vsearchPreResults, vsearchPostResults, "Vsearch", "cyan4")

plotFun(pathwayPreResults, pathwayPostResults, "Pathway", "khaki")

dev.off()

