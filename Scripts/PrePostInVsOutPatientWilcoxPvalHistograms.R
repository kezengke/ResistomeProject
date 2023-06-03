#in vs out patient wilcoxon test for pre and post patient histograms.
rm(list = ls())

#wilcox results
brackenPreResults<-read.csv("Wilcox/brackenPreInOut.csv", row.names = 1)
brackenPostResults<-read.csv("Wilcox/brackenPostInOut.csv", row.names = 1)
amrPreResults<-read.csv("Wilcox/amrPreInOut.csv", row.names = 1)
amrPostResults<-read.csv("Wilcox/amrPostInOut.csv", row.names = 1)
rgiPreResults<-read.csv("Wilcox/rgiPreInOut.csv", row.names = 1)
rgiPostResults<-read.csv("Wilcox/rgiPostInOut.csv", row.names = 1)
vsearchPreResults<-read.csv("Wilcox/vsearchPreInOut.csv", row.names = 1)
vsearchPostResults<-read.csv("Wilcox/vsearchPostInOut.csv", row.names = 1)

pdf("Plots/PrePostInVsOutPatientWilcoxonPvalHistograms.pdf", width=12, height=24)
par(mfrow=c(4, 2))
par(mar=c(5, 6, 4, 1)+.1)
#bracken
hist(brackenPreResults$brackenPreWilcox_p, breaks=seq(0, 1, 0.05), xlab = "p-value", col = "tan2",
     main = "Bracken(Species) Pre In vs Out Wilcoxon P-vals", 
     xlim=c(0,1), cex.lab = 1.5, cex.main = 1.7, cex.axis = 1.4)
abline(v=0.05, col=gray(.5), lty=2)
mtext(at=0.05, side=3, text=0.05, col=gray(.5))

hist(brackenPostResults$brackenPostWilcox_p, breaks=seq(0, 1, 0.05), xlab = "p-value", col = "tan2",
     main = "Bracken(Species) Post In vs Out Wilcoxon P-vals", 
     xlim=c(0,1), cex.lab = 1.5, cex.main = 1.7, cex.axis = 1.4)
abline(v=0.05, col=gray(.5), lty=2)
mtext(at=0.05, side=3, text=0.05, col=gray(.5))
#amr
hist(amrPreResults$amrPreWilcox_p, breaks=seq(0, 1, 0.05), xlab = "p-value", col = "coral3",
     main = "AMR Pre In vs Out Wilcoxon P-vals", 
     xlim=c(0,1), cex.lab = 1.5, cex.main = 1.7, cex.axis = 1.4)
abline(v=0.05, col=gray(.5), lty=2)
mtext(at=0.05, side=3, text=0.05, col=gray(.5))

hist(amrPostResults$amrPostWilcox_p, breaks=seq(0, 1, 0.05), xlab = "p-value", col = "coral3",
     main = "AMR Post In vs Out Wilcoxon P-vals", 
     xlim=c(0,1), cex.lab = 1.5, cex.main = 1.7, cex.axis = 1.4)
abline(v=0.05, col=gray(.5), lty=2)
mtext(at=0.05, side=3, text=0.05, col=gray(.5))
#rgi
hist(rgiPreResults$rgiPreWilcox_p, breaks=seq(0, 1, 0.05), xlab = "p-value", col = "cornflowerblue",
     main = "RGI Pre In vs Out Wilcoxon P-vals", 
     xlim=c(0,1), cex.lab = 1.5, cex.main = 1.7, cex.axis = 1.4)
abline(v=0.05, col=gray(.5), lty=2)
mtext(at=0.05, side=3, text=0.05, col=gray(.5))

hist(rgiPostResults$rgiPostWilcox_p, breaks=seq(0, 1, 0.05), xlab = "p-value", col = "cornflowerblue",
     main = "RGI Post In vs Out Wilcoxon P-vals", 
     xlim=c(0,1), cex.lab = 1.5, cex.main = 1.7, cex.axis = 1.4)
abline(v=0.05, col=gray(.5), lty=2)
mtext(at=0.05, side=3, text=0.05, col=gray(.5))
#vsearch
hist(vsearchPreResults$vsearchPreWilcox_p, breaks=seq(0, 1, 0.05), xlab = "p-value", col = "olivedrab4",
     main = "vsearch Pre In vs Out Wilcoxon P-vals", 
     xlim=c(0,1), cex.lab = 1.5, cex.main = 1.7, cex.axis = 1.4)
abline(v=0.05, col=gray(.5), lty=2)
mtext(at=0.05, side=3, text=0.05, col=gray(.5))

hist(vsearchPostResults$vsearchPostWilcox_p, breaks=seq(0, 1, 0.05), xlab = "p-value", col = "olivedrab4",
     main = "vsearch Post In vs Out Wilcoxon P-vals", 
     xlim=c(0,1), cex.lab = 1.5, cex.main = 1.7, cex.axis = 1.4)
abline(v=0.05, col=gray(.5), lty=2)
mtext(at=0.05, side=3, text=0.05, col=gray(.5))

dev.off()