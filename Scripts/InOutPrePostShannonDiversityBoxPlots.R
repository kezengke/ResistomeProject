#in vs out patient pre post shannon diversity box plots
rm(list = ls())
library("vegan")

metaData<-read.csv("metaWithBins.csv", header = T, row.names = 1)

#counts tables
brackenT<-read.csv("CountsTables/brackenProcessed.csv", header = T, row.names = 1, check.names = F)
amrT<-read.csv("CountsTables/amrProcessed.csv", header = T, row.names = 1, check.names = F)
rgiT<-read.csv("CountsTables/rgiProcessed.csv", header = T, row.names = 1, check.names = F)
vsearchT<-read.csv("CountsTables/vsearchProcessed.csv", header = T, row.names = 1, check.names = F)
genusT<-read.csv("CountsTables/genusProcessed.csv", header = T, row.names = 1, check.names = F)

patients<-unique(metaData$ID)

pres<-vector()
posts<-vector()
for (i in 1:length(patients)) {
  patient<-metaData[metaData$ID == patients[i], , drop = F]
  pres[i]<-rownames(patient)[patient$bins == "PRE"]
  posts[i]<-rownames(patient)[which.min(abs(patient$Timepoint - 28))]
}

brackenPre<-brackenT[, pres, drop = F]
brackenPost<-brackenT[, posts, drop = F]
amrPre<-amrT[, pres, drop = F]
amrPost<-amrT[, intersect(colnames(amrT), posts), drop = F]
rgiPre<-rgiT[, pres, drop = F]
rgiPost<-rgiT[, posts, drop = F]
vsearchPre<-vsearchT[, pres, drop = F]
vsearchPost<-vsearchT[, posts, drop = F]
genusPre<-genusT[, pres, drop = F]
genusPost<-genusT[, posts, drop = F]

metaBRACKENpre<-metaData[colnames(brackenPre), , drop = F]
metaBRACKENpost<-metaData[colnames(brackenPost), , drop = F]
metaAMRpre<-metaData[colnames(amrPre), , drop = F]
metaAMRpost<-metaData[colnames(amrPost), , drop = F]
metaRGIpre<-metaData[colnames(rgiPre), , drop = F]
metaRGIpost<-metaData[colnames(rgiPost), , drop = F]
metaVSEARCHpre<-metaData[colnames(vsearchPre), , drop = F]
metaVSEARCHpost<-metaData[colnames(vsearchPost), , drop = F]
metaGENUSpre<-metaData[colnames(genusPre), , drop = F]
metaGENUSpost<-metaData[colnames(genusPost), , drop = F]

pdf("Plots/ShannonDiversityBoxPlotsForInVsOutPatientsInPreAndD28.pdf", width=7, height=10.5)
par(mfrow=c(3, 2))
par(mar=c(5,6,4,1)+.1)
#bracken
T1pre<-brackenPre[, which(metaBRACKENpre$ptInOut == "Inpatient"), drop = F]
T2pre<-brackenPre[, which(metaBRACKENpre$ptInOut == "Outpatient"), drop = F]

H1pre<-diversity(T1pre, MARGIN = 2)
H2pre<-diversity(T2pre, MARGIN = 2)
H3pre<-diversity(brackenPre, MARGIN = 2)

prewilcoxP<-wilcox.test(H3pre~metaBRACKENpre$ptInOut)$p.value
x<-list("Inpatient"=H1pre, "Outpatient"=H2pre)
boxplot(x, outline = F,
        cex.lab = 1.5, cex.axis = 1.4, cex.main = 1.5, 
        col=c("olivedrab3", "cornflowerblue"),
        ylab = "Shannon Diversity",
        xlab = "PRE",
        main = paste0("(Species)P-value=",signif(prewilcoxP, 4)))
stripchart(x, method="jitter", 
           vertical=T, pch=19, add=T)

T1post<-brackenPost[, which(metaBRACKENpost$ptInOut == "Inpatient"), drop = F]
T2post<-brackenPost[, which(metaBRACKENpost$ptInOut == "Outpatient"), drop = F]

H1post<-diversity(T1post, MARGIN = 2)
H2post<-diversity(T2post, MARGIN = 2)
H3post<-diversity(brackenPost, MARGIN = 2)

postwilcoxP<-wilcox.test(H3post~metaBRACKENpost$ptInOut)$p.value
x<-list("Inpatient"=H1post, "Outpatient"=H2post)
boxplot(x, outline = F,
        cex.lab = 1.5, cex.axis = 1.4, cex.main = 1.5, 
        col=c("olivedrab3", "cornflowerblue"),
        ylab = "Shannon Diversity",
        xlab = "D28",
        main = paste0("(Species)P-value=",signif(postwilcoxP, 4)))
stripchart(x, method="jitter", 
           vertical=T, pch=19, add=T)

#genus
T1pre<-genusPre[, which(metaGENUSpre$ptInOut == "Inpatient"), drop = F]
T2pre<-genusPre[, which(metaGENUSpre$ptInOut == "Outpatient"), drop = F]

H1pre<-diversity(T1pre, MARGIN = 2)
H2pre<-diversity(T2pre, MARGIN = 2)
H3pre<-diversity(genusPre, MARGIN = 2)

prewilcoxP<-wilcox.test(H3pre~metaGENUSpre$ptInOut)$p.value
x<-list("Inpatient"=H1pre, "Outpatient"=H2pre)
boxplot(x, outline = F,
        cex.lab = 1.5, cex.axis = 1.4, cex.main = 1.5, 
        col=c("olivedrab3", "cornflowerblue"),
        ylab = "Shannon Diversity",
        xlab = "PRE",
        main = paste0("(Genus)P-value=",signif(prewilcoxP, 4)))
stripchart(x, method="jitter", 
           vertical=T, pch=19, add=T)

T1post<-genusPost[, which(metaGENUSpost$ptInOut == "Inpatient"), drop = F]
T2post<-genusPost[, which(metaGENUSpost$ptInOut == "Outpatient"), drop = F]

H1post<-diversity(T1post, MARGIN = 2)
H2post<-diversity(T2post, MARGIN = 2)
H3post<-diversity(genusPost, MARGIN = 2)

postwilcoxP<-wilcox.test(H3post~metaGENUSpost$ptInOut)$p.value
x<-list("Inpatient"=H1post, "Outpatient"=H2post)
boxplot(x, outline = F,
        cex.lab = 1.5, cex.axis = 1.4, cex.main = 1.5, 
        col=c("olivedrab3", "cornflowerblue"),
        ylab = "Shannon Diversity",
        xlab = "D28",
        main = paste0("(Genus)P-value=",signif(postwilcoxP, 4)))
stripchart(x, method="jitter", 
           vertical=T, pch=19, add=T)

#amr
T1pre<-amrPre[, which(metaAMRpre$ptInOut == "Inpatient"), drop = F]
T2pre<-amrPre[, which(metaAMRpre$ptInOut == "Outpatient"), drop = F]

H1pre<-diversity(T1pre, MARGIN = 2)
H2pre<-diversity(T2pre, MARGIN = 2)
H3pre<-diversity(amrPre, MARGIN = 2)

prewilcoxP<-wilcox.test(H3pre~metaAMRpre$ptInOut)$p.value
x<-list("Inpatient"=H1pre, "Outpatient"=H2pre)
boxplot(x, outline = F,
        cex.lab = 1.5, cex.axis = 1.4, cex.main = 1.5, 
        col=c("olivedrab3", "cornflowerblue"),
        ylab = "Shannon Diversity",
        xlab = "PRE",
        main = paste0("(AMR)P-value=",signif(prewilcoxP, 4)))
stripchart(x, method="jitter", 
           vertical=T, pch=19, add=T)

T1post<-amrPost[, which(metaAMRpost$ptInOut == "Inpatient"), drop = F]
T2post<-amrPost[, which(metaAMRpost$ptInOut == "Outpatient"), drop = F]

H1post<-diversity(T1post, MARGIN = 2)
H2post<-diversity(T2post, MARGIN = 2)
H3post<-diversity(amrPost, MARGIN = 2)

postwilcoxP<-wilcox.test(H3post~metaAMRpost$ptInOut)$p.value
x<-list("Inpatient"=H1post, "Outpatient"=H2post)
boxplot(x, outline = F,
        cex.lab = 1.5, cex.axis = 1.4, cex.main = 1.5, 
        col=c("olivedrab3", "cornflowerblue"),
        ylab = "Shannon Diversity",
        xlab = "D28",
        main = paste0("(AMR)P-value=",signif(postwilcoxP, 4)))
stripchart(x, method="jitter", 
           vertical=T, pch=19, add=T)

#rgi
T1pre<-rgiPre[, which(metaRGIpre$ptInOut == "Inpatient"), drop = F]
T2pre<-rgiPre[, which(metaRGIpre$ptInOut == "Outpatient"), drop = F]

H1pre<-diversity(T1pre, MARGIN = 2)
H2pre<-diversity(T2pre, MARGIN = 2)
H3pre<-diversity(rgiPre, MARGIN = 2)

prewilcoxP<-wilcox.test(H3pre~metaRGIpre$ptInOut)$p.value
x<-list("Inpatient"=H1pre, "Outpatient"=H2pre)
boxplot(x, outline = F,
        cex.lab = 1.5, cex.axis = 1.4, cex.main = 1.5, 
        col=c("olivedrab3", "cornflowerblue"),
        ylab = "Shannon Diversity",
        xlab = "PRE",
        main = paste0("(RGI)P-value=",signif(prewilcoxP, 4)))
stripchart(x, method="jitter", 
           vertical=T, pch=19, add=T)

T1post<-rgiPost[, which(metaRGIpost$ptInOut == "Inpatient"), drop = F]
T2post<-rgiPost[, which(metaRGIpost$ptInOut == "Outpatient"), drop = F]

H1post<-diversity(T1post, MARGIN = 2)
H2post<-diversity(T2post, MARGIN = 2)
H3post<-diversity(rgiPost, MARGIN = 2)

postwilcoxP<-wilcox.test(H3post~metaRGIpost$ptInOut)$p.value
x<-list("Inpatient"=H1post, "Outpatient"=H2post)
boxplot(x, outline = F,
        cex.lab = 1.5, cex.axis = 1.4, cex.main = 1.5, 
        col=c("olivedrab3", "cornflowerblue"),
        ylab = "Shannon Diversity",
        xlab = "D28",
        main = paste0("(RGI)P-value=",signif(postwilcoxP, 4)))
stripchart(x, method="jitter", 
           vertical=T, pch=19, add=T)

#vsearch
T1pre<-vsearchPre[, which(metaVSEARCHpre$ptInOut == "Inpatient"), drop = F]
T2pre<-vsearchPre[, which(metaVSEARCHpre$ptInOut == "Outpatient"), drop = F]

H1pre<-diversity(T1pre, MARGIN = 2)
H2pre<-diversity(T2pre, MARGIN = 2)
H3pre<-diversity(vsearchPre, MARGIN = 2)

prewilcoxP<-wilcox.test(H3pre~metaVSEARCHpre$ptInOut)$p.value
x<-list("Inpatient"=H1pre, "Outpatient"=H2pre)
boxplot(x, outline = F,
        cex.lab = 1.5, cex.axis = 1.4, cex.main = 1.5, 
        col=c("olivedrab3", "cornflowerblue"),
        ylab = "Shannon Diversity",
        xlab = "PRE",
        main = paste0("(vsearch)P-value=",signif(prewilcoxP, 4)))
stripchart(x, method="jitter", 
           vertical=T, pch=19, add=T)

T1post<-vsearchPost[, which(metaVSEARCHpost$ptInOut == "Inpatient"), drop = F]
T2post<-vsearchPost[, which(metaVSEARCHpost$ptInOut == "Outpatient"), drop = F]

H1post<-diversity(T1post, MARGIN = 2)
H2post<-diversity(T2post, MARGIN = 2)
H3post<-diversity(vsearchPost, MARGIN = 2)

postwilcoxP<-wilcox.test(H3post~metaVSEARCHpost$ptInOut)$p.value
x<-list("Inpatient"=H1post, "Outpatient"=H2post)
boxplot(x, outline = F,
        cex.lab = 1.5, cex.axis = 1.4, cex.main = 1.5, 
        col=c("olivedrab3", "cornflowerblue"),
        ylab = "Shannon Diversity",
        xlab = "D28",
        main = paste0("(vsearch)P-value=",signif(postwilcoxP, 4)))
stripchart(x, method="jitter", 
           vertical=T, pch=19, add=T)

dev.off()

