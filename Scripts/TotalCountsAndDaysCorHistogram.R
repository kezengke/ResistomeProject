rm(list = ls())
library(stringr)
pdf("Plots/CorrelationHistograms.pdf", width=10, height=15)
par(mfrow=c(3,2))

makeMeta<-read.csv("Tessa_duke_newPipeline_count_table[1].csv", sep = ",", header = T, row.names = 1)
rownames(makeMeta)[grepl("D20726D62", rownames(makeMeta))]
RownameToRemove<-c("D20726D62_CGAGGCTG-ATTAGACG_S91_L003", "D21685NPE_TAAGGCGA-GCGTAAGA_S249_L004")
makeMeta<-makeMeta[!(rownames(makeMeta) %in% RownameToRemove), ]
rownames(makeMeta)<-sapply(str_split(rownames(makeMeta), "_",  n = 4), `[`, 1)#get rid of random -s
rownames(makeMeta)<-gsub("-", "", rownames(makeMeta))

sampleNames<-rownames(makeMeta)

#AMR
makeAMR<-read.delim("CountsTables/AMR_counts.tsv", sep = "\t", header = T, row.names = 1)
colnames(makeAMR)<-gsub("pre", "PRE", colnames(makeAMR), ignore.case = T) #switch all pre to cap
colnames(makeAMR)<-sapply(str_split(colnames(makeAMR), "_", n = 4 ), `[`, 1) #get rid of useless info
colnames(makeAMR)<-gsub("\\.", "", colnames(makeAMR))#get rid of random .s

makeAMR<-makeAMR[, intersect(sampleNames, colnames(makeAMR))]
ID<-unique(paste0("D", sapply(str_split(colnames(makeAMR)[!grepl("PRE", colnames(makeAMR))], "D", n = 3 ), `[`, 2)))

pearsonCor<-c()
spearmanRho<-c()

for (i in 1:length(ID)) {
  allNames<-colnames(makeAMR)[grepl(ID[i], colnames(makeAMR))]
  days<-sort(as.numeric(sapply(str_split(colnames(makeAMR)[grepl(ID[i], colnames(makeAMR))], "D", n = 3), `[`, 3)))
  orderedNames<-paste(ID[i], days, sep = "D")
  if (sum(grepl("PRE", allNames, ignore.case = T))>0){
    orderedNames<-append(paste0(ID[i], "PRE"), orderedNames)
    days<-append(0, days)
  }
  
  if(length(days)<3)
    next
  pearsonCor<-c(pearsonCor, cor.test(days, colSums(makeAMR[, orderedNames]), method = 'pearson')$estimate)
  spearmanRho<-c(spearmanRho, cor.test(days, colSums(makeAMR[, orderedNames]), method = 'spearman')$estimate)
}

hist(pearsonCor, breaks=seq(-1, 1, 0.1), col = "coral3", 
     xlab = "Pearson's Cor", main="AMR", 
    cex.lab = 1.5, cex.main = 1.7, cex.axis = 1.4)
hist(spearmanRho, breaks=seq(-1, 1, 0.1), col = "coral3", 
     xlab = "Spearman's Rho", main="AMR", 
     cex.lab = 1.5, cex.main = 1.7, cex.axis = 1.4)

#RGI
makeRGI<-read.delim("CountsTables/RGI_counts.tsv", sep = "\t", header = T, row.names = 1)
colnames(makeRGI)<-gsub("pre", "PRE", colnames(makeRGI), ignore.case = T) #switch all pre to cap
colnames(makeRGI)<-sapply(str_split(colnames(makeRGI), "_", n = 4 ), `[`, 1) #get rid of useless info
colnames(makeRGI)<-gsub("\\.", "", colnames(makeRGI))#get rid of random .s

makeRGI<-makeRGI[, intersect(sampleNames, colnames(makeRGI))]
ID<-unique(paste0("D", sapply(str_split(colnames(makeRGI)[!grepl("PRE", colnames(makeRGI))], "D", n = 3 ), `[`, 2)))

pearsonCor<-c()
spearmanRho<-c()

for (i in 1:length(ID)) {
  allNames<-colnames(makeRGI)[grepl(ID[i], colnames(makeRGI))]
  days<-sort(as.numeric(sapply(str_split(colnames(makeRGI)[grepl(ID[i], colnames(makeRGI))], "D", n = 3), `[`, 3)))
  orderedNames<-paste(ID[i], days, sep = "D")
  if (sum(grepl("PRE", allNames, ignore.case = T))>0){
    orderedNames<-append(paste0(ID[i], "PRE"), orderedNames)
    days<-append(0, days)
  }
  
  if(length(days)<3)
    next
  pearsonCor<-c(pearsonCor, cor.test(days, colSums(makeRGI[, orderedNames]), method = 'pearson')$estimate)
  spearmanRho<-c(spearmanRho, cor.test(days, colSums(makeRGI[, orderedNames]), method = 'spearman')$estimate)
}

hist(pearsonCor, breaks=seq(-1, 1, 0.1), col = "cornflowerblue", 
     xlab = "Pearson's Cor", main="RGI", 
     cex.lab = 1.5, cex.main = 1.7, cex.axis = 1.4)
hist(spearmanRho, breaks=seq(-1, 1, 0.1), col = "cornflowerblue", 
     xlab = "Spearman's Rho", main="RGI", 
     cex.lab = 1.5, cex.main = 1.7, cex.axis = 1.4)

#Vsearch
makeVsearch<-read.delim("CountsTables/vsearch_counts.tsv", sep = "\t", header = T, row.names = 1)
colnames(makeVsearch)<-gsub("pre", "PRE", colnames(makeVsearch), ignore.case = T) #switch all pre to cap
colnames(makeVsearch)<-sapply(str_split(colnames(makeVsearch), "_", n = 4 ), `[`, 1) #get rid of useless info
colnames(makeVsearch)<-gsub("\\.", "", colnames(makeVsearch))#get rid of random .s

makeVsearch<-makeVsearch[, intersect(sampleNames, colnames(makeVsearch))]
ID<-unique(paste0("D", sapply(str_split(colnames(makeVsearch)[!grepl("PRE", colnames(makeVsearch))], "D", n = 3 ), `[`, 2)))

pearsonCor<-c()
spearmanRho<-c()

for (i in 1:length(ID)) {
  allNames<-colnames(makeVsearch)[grepl(ID[i], colnames(makeVsearch))]
  days<-sort(as.numeric(sapply(str_split(colnames(makeVsearch)[grepl(ID[i], colnames(makeVsearch))], "D", n = 3), `[`, 3)))
  orderedNames<-paste(ID[i], days, sep = "D")
  if (sum(grepl("PRE", allNames, ignore.case = T))>0){
    orderedNames<-append(paste0(ID[i], "PRE"), orderedNames)
    days<-append(0, days)
  }
  
  if(length(days)<3)
    next
  pearsonCor<-c(pearsonCor, cor.test(days, colSums(makeVsearch[, orderedNames]), method = 'pearson')$estimate)
  spearmanRho<-c(spearmanRho, cor.test(days, colSums(makeVsearch[, orderedNames]), method = 'spearman')$estimate)
}

hist(pearsonCor, breaks=seq(-1, 1, 0.1), col = "olivedrab4", 
     xlab = "Pearson's Cor", main="Vsearch", 
     cex.lab = 1.5, cex.main = 1.7, cex.axis = 1.4)
hist(spearmanRho, breaks=seq(-1, 1, 0.1), col = "olivedrab4", 
     xlab = "Spearman's Rho", main="Vsearch", 
     cex.lab = 1.5, cex.main = 1.7, cex.axis = 1.4)

dev.off()
