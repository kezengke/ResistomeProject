rm(list = ls())
library(stringr)
library(ggplot2)
makeMeta<-read.csv("Tessa_duke_newPipeline_count_table[1].csv", sep = ",", header = T, row.names = 1)
rownames(makeMeta)[grepl("D20726D62", rownames(makeMeta))]
RownameToRemove<-c("D20726D62_CGAGGCTG-ATTAGACG_S91_L003", "D21685NPE_TAAGGCGA-GCGTAAGA_S249_L004")
makeMeta<-makeMeta[!(rownames(makeMeta) %in% RownameToRemove), ]
rownames(makeMeta)<-sapply(str_split(rownames(makeMeta), "_",  n = 4), `[`, 1)#get rid of random -s
rownames(makeMeta)<-gsub("-", "", rownames(makeMeta))

sampleNames<-rownames(makeMeta)

makeAMR<-read.delim("CountsTables/AMR_counts.tsv", sep = "\t", header = T, row.names = 1)
colnames(makeAMR)<-gsub("pre", "PRE", colnames(makeAMR), ignore.case = T) #switch all pre to cap
colnames(makeAMR)<-sapply(str_split(colnames(makeAMR), "_", n = 4 ), `[`, 1) #get rid of useless info
colnames(makeAMR)<-gsub("\\.", "", colnames(makeAMR))#get rid of random .s

makeAMR<-makeAMR[, intersect(sampleNames, colnames(makeAMR))]
ID<-unique(paste0("D", sapply(str_split(colnames(makeAMR)[!grepl("PRE", colnames(makeAMR))], "D", n = 3 ), `[`, 2)))

allNames<-colnames(makeAMR)[grepl(ID[1], colnames(makeAMR))]
days<-sort(as.numeric(sapply(str_split(colnames(makeAMR)[grepl(ID[1], colnames(makeAMR))], "D", n = 3), `[`, 3)))
orderedNames<-paste(ID[1], days, sep = "D")
if (sum(grepl("PRE", allNames, ignore.case = T))>0)
  orderedNames<-append(paste0(ID[1], "PRE"), orderedNames)
View(makeAMR[, orderedNames])


ggplot(makeAMR[, orderedNames], aes(fill=rownames(makeAMR[, orderedNames]), y=makeAMR[, orderedNames], x=colnames(makeAMR[, orderedNames]))) + 
  geom_bar(position="stack", stat="identity")

for (i in 1:nrow(makeAMR)) {
  plot(unlist(makeAMR[i, orderedNames]), col = "coral3", pch = 19, xaxt = "n", 
       xlab = "Days", ylab = "Counts", 
       main = paste(ID[1], ":", rownames(makeAMR[i, orderedNames])))
  axis(1, at = 1:length(orderedNames), labels = append("PRE", paste0("D", days)),  cex = 0.3, adj = 1)
}

makeRGI<-read.delim("CountsTables/RGI_counts.tsv", sep = "\t", header = T, row.names = 1)
# sapply(str_split(colnames(makeRGI), "_",  n = 4), `[`, 1)
# length(unique(sapply(str_split(colnames(makeRGI), "_",  n = 4), `[`, 1)))
colnames(makeRGI)<-sapply(str_split(colnames(makeRGI), "_", n = 4 ), `[`, 1)

makeVSEARCH<-read.delim("CountsTables/vsearch_counts.tsv", sep = "\t", header = T, row.names = 1)
# sapply(str_split(colnames(makeRGI), "_",  n = 4), `[`, 1)
# length(unique(sapply(str_split(colnames(makeRGI), "_",  n = 4), `[`, 1)))
colnames(makeVSEARCH)<-sapply(str_split(colnames(makeVSEARCH), "_", n = 4 ), `[`, 1)

