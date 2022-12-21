rm(list = ls())
library(stringr)
# makeMeta<-read.csv("Tessa_duke_newPipeline_count_table[1].csv", sep = ",", header = T, row.names = 1)
# rownames(makeMeta)<-sapply(str_split(rownames(makeMeta), "_",  n = 4), `[`, 3)

makeAMR<-read.delim("CountsTables/AMR_counts.tsv", sep = "\t", header = T, row.names = 1)
# sapply(str_split(colnames(makeAMR), "_",  n = 4), `[`, 1)
# length(unique(sapply(str_split(colnames(makeAMR), "_",  n = 4), `[`, 1)))
colnames(makeAMR)<-str_replace(colnames(makeAMR), regex("pre", ignore_case = T), "PRE") #switch all pre to cap
colnames(makeAMR)<-sapply(str_split(colnames(makeAMR), "_", n = 4 ), `[`, 1) #get rid off useless info

sampleNames<-unique(sapply(str_split(colnames(makeAMR), ".D", n = 2 ), `[`, 1))
# sampleNames<-unique(sapply(str_split(sapply(str_split(colnames(makeAMR), ".D", n = 2 ), `[`, 1), regex(".PRE", ignore_case = T), n = 2 ), `[`, 1))
# for (i in 1:length(sampleNames)) {
#   print(colnames(makeAMR)[grepl(sampleNames[i], colnames(makeAMR))])
# }

allNames<-colnames(makeAMR)[grepl(sampleNames[1], colnames(makeAMR))]
days<-sort(as.numeric(sapply(str_split(colnames(makeAMR)[grepl(sampleNames[1], colnames(makeAMR))], ".D", n = 2), `[`, 2)))
# allDays<-sapply(str_split(colnames(makeAMR)[grepl(sampleNames[1], colnames(makeAMR))], ".D", n = 2), `[`, 2)
orderedNames<-paste0(sampleNames[1], days, sep = "D")
if (sum(grepl("PRE", allNames, ignore.case = T))>0)
  days<-append(paste0(sampleNames[1], "PRE"), days)



makeRGI<-read.delim("CountsTables/RGI_counts.tsv", sep = "\t", header = T, row.names = 1)
# sapply(str_split(colnames(makeRGI), "_",  n = 4), `[`, 1)
# length(unique(sapply(str_split(colnames(makeRGI), "_",  n = 4), `[`, 1)))
colnames(makeRGI)<-sapply(str_split(colnames(makeRGI), "_", n = 4 ), `[`, 1)

makeVSEARCH<-read.delim("CountsTables/vsearch_counts.tsv", sep = "\t", header = T, row.names = 1)
# sapply(str_split(colnames(makeRGI), "_",  n = 4), `[`, 1)
# length(unique(sapply(str_split(colnames(makeRGI), "_",  n = 4), `[`, 1)))
colnames(makeVSEARCH)<-sapply(str_split(colnames(makeVSEARCH), "_", n = 4 ), `[`, 1)

