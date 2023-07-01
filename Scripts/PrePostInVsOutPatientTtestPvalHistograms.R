#in vs out patient t-test for pre and post patient histograms.
rm(list = ls())

#wilcox results
brackenPreResults<-read.csv("Ttest/brackenPreInOut.csv", row.names = 1)
brackenPostResults<-read.csv("Ttest/brackenPostInOut.csv", row.names = 1)
amrPreResults<-read.csv("Ttest/amrPreInOut.csv", row.names = 1)
amrPostResults<-read.csv("Ttest/amrPostInOut.csv", row.names = 1)
rgiPreResults<-read.csv("Ttest/rgiPreInOut.csv", row.names = 1)
rgiPostResults<-read.csv("Ttest/rgiPostInOut.csv", row.names = 1)
vsearchPreResults<-read.csv("Ttest/vsearchPreInOut.csv", row.names = 1)
vsearchPostResults<-read.csv("Ttest/vsearchPostInOut.csv", row.names = 1)
genusPreResults<-read.csv("Ttest/genusPreInOut.csv", row.names = 1)
genusPostResults<-read.csv("Ttest/genusPostInOut.csv", row.names = 1)

pdf("Plots/PrePostInVsOutPatientTtestPvalHistograms.pdf", width=12, height=12)
par(mfrow=c(2, 2))
par(mar=c(5, 6, 4, 1)+.1)
#bracken
hist(brackenPreResults$brackenPreTtest_p, breaks=seq(0, 1, 0.05), xlab = "p-value", col = "tan2",
     main = "Bracken(Species) Pre In vs Out Ttest P-vals", 
     xlim=c(0,1), cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.4)
abline(v=0.05, col=gray(.5), lty=2)
mtext(at=0.05, side=3, text=0.05, col=gray(.5))

hist(brackenPostResults$brackenPostTtest_p, breaks=seq(0, 1, 0.05), xlab = "p-value", col = "tan2",
     main = "Bracken(Species) Post In vs Out Ttest P-vals", 
     xlim=c(0,1), cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.4)
abline(v=0.05, col=gray(.5), lty=2)
mtext(at=0.05, side=3, text=0.05, col=gray(.5))

#genus
hist(genusPreResults$genusPreTtest_p, breaks=seq(0, 1, 0.05), xlab = "p-value", col = "orchid3",
     main = "Bracken(Genus) Pre In vs Out Ttest P-vals", 
     xlim=c(0,1), cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.4)
abline(v=0.05, col=gray(.5), lty=2)
mtext(at=0.05, side=3, text=0.05, col=gray(.5))

hist(genusPostResults$genusPostTtest_p, breaks=seq(0, 1, 0.05), xlab = "p-value", col = "orchid3",
     main = "Bracken(Genus) Post In vs Out Ttest P-vals", 
     xlim=c(0,1), cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.4)
abline(v=0.05, col=gray(.5), lty=2)
mtext(at=0.05, side=3, text=0.05, col=gray(.5))

#amr
hist(amrPreResults$amrPreTtest_p, breaks=seq(0, 1, 0.05), xlab = "p-value", col = "coral3",
     main = "AMR Pre In vs Out Ttest P-vals", 
     xlim=c(0,1), cex.lab = 1.5, cex.main = 1.7, cex.axis = 1.4)
abline(v=0.05, col=gray(.5), lty=2)
mtext(at=0.05, side=3, text=0.05, col=gray(.5))

hist(amrPostResults$amrPostTtest_p, breaks=seq(0, 1, 0.05), xlab = "p-value", col = "coral3",
     main = "AMR Post In vs Out Ttest P-vals", 
     xlim=c(0,1), cex.lab = 1.5, cex.main = 1.7, cex.axis = 1.4)
abline(v=0.05, col=gray(.5), lty=2)
mtext(at=0.05, side=3, text=0.05, col=gray(.5))
#rgi
hist(rgiPreResults$rgiPreTtest_p, breaks=seq(0, 1, 0.05), xlab = "p-value", col = "cornflowerblue",
     main = "RGI Pre In vs Out Ttest P-vals", 
     xlim=c(0,1), cex.lab = 1.5, cex.main = 1.7, cex.axis = 1.4)
abline(v=0.05, col=gray(.5), lty=2)
mtext(at=0.05, side=3, text=0.05, col=gray(.5))

hist(rgiPostResults$rgiPostTtest_p, breaks=seq(0, 1, 0.05), xlab = "p-value", col = "cornflowerblue",
     main = "RGI Post In vs Out Ttest P-vals", 
     xlim=c(0,1), cex.lab = 1.5, cex.main = 1.7, cex.axis = 1.4)
abline(v=0.05, col=gray(.5), lty=2)
mtext(at=0.05, side=3, text=0.05, col=gray(.5))
#vsearch
hist(vsearchPreResults$vsearchPreTtest_p, breaks=seq(0, 1, 0.05), xlab = "p-value", col = "olivedrab4",
     main = "vsearch Pre In vs Out Ttest P-vals", 
     xlim=c(0,1), cex.lab = 1.5, cex.main = 1.7, cex.axis = 1.4)
abline(v=0.05, col=gray(.5), lty=2)
mtext(at=0.05, side=3, text=0.05, col=gray(.5))

hist(vsearchPostResults$vsearchPostTtest_p, breaks=seq(0, 1, 0.05), xlab = "p-value", col = "olivedrab4",
     main = "vsearch Post In vs Out Ttest P-vals", 
     xlim=c(0,1), cex.lab = 1.5, cex.main = 1.7, cex.axis = 1.4)
abline(v=0.05, col=gray(.5), lty=2)
mtext(at=0.05, side=3, text=0.05, col=gray(.5))

dev.off()