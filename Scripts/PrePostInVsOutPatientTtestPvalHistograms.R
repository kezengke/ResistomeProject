#in vs out patient t-test for pre and post patient histograms.
rm(list = ls())

plotFun <- function(resPre, resPost, resType, color) {
        hist(resPre[, 1], breaks=seq(0, 1, 0.05), xlab = "p-value", col = color,
             main = paste0("(", resType, ")", "Pre In vs Out Ttest P-vals"), 
             xlim=c(0,1), cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.4)
        abline(v=0.05, col=gray(.5), lty=2)
        mtext(at=0.05, side=3, text=0.05, col=gray(.5))
        
        hist(resPost[, 1], breaks=seq(0, 1, 0.05), xlab = "p-value", col = color,
             main = paste0("(", resType, ")", "Post In vs Out Ttest P-vals"), 
             xlim=c(0,1), cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.4)
        abline(v=0.05, col=gray(.5), lty=2)
        mtext(at=0.05, side=3, text=0.05, col=gray(.5))
}

#wilcox results
speciesPreResults<-read.csv("Ttest/speciesPreInOut.csv", row.names = 1)
speciesPostResults<-read.csv("Ttest/speciesPostInOut.csv", row.names = 1)
genusPreResults<-read.csv("Ttest/genusPreInOut.csv", row.names = 1)
genusPostResults<-read.csv("Ttest/genusPostInOut.csv", row.names = 1)
phylumPreResults<-read.csv("Ttest/phylumPreInOut.csv", row.names = 1)
phylumPostResults<-read.csv("Ttest/phylumPostInOut.csv", row.names = 1)

amrPreResults<-read.csv("Ttest/amrPreInOut.csv", row.names = 1)
amrPostResults<-read.csv("Ttest/amrPostInOut.csv", row.names = 1)
rgiPreResults<-read.csv("Ttest/rgiPreInOut.csv", row.names = 1)
rgiPostResults<-read.csv("Ttest/rgiPostInOut.csv", row.names = 1)
vsearchPreResults<-read.csv("Ttest/vsearchPreInOut.csv", row.names = 1)
vsearchPostResults<-read.csv("Ttest/vsearchPostInOut.csv", row.names = 1)

pathwayPreResults<-read.csv("Ttest/pathwayPreInOut.csv", row.names = 1)
pathwayPostResults<-read.csv("Ttest/pathwayPostInOut.csv", row.names = 1)

pdf("Plots/PrePostInVsOutPatientTtestPvalHistograms.pdf", width=12, height=18)
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