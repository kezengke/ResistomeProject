rm(list = ls())
dukeSamples<-read.csv("dukeSamples.csv", header = T, sep = ",")
dukeSamples$Sample<-gsub("-", "", dukeSamples$Sample)
dukeSamples$Sample<-gsub(" ", "", dukeSamples$Sample)

#AMR
amrT<-read.delim("CountsTables/AMR_counts.tsv", sep = "\t", header = T, row.names = 1)
colnames(amrT)<-gsub("pre", "PRE", colnames(amrT), ignore.case = T) #switch all pre to cap
colnames(amrT)<-sapply(str_split(colnames(amrT), "_", n = 4 ), `[`, 1) #get rid of useless info
colnames(amrT)<-gsub("\\.", "", colnames(amrT))#get rid of random .s

length(intersect(dukeSamples$Sample, colnames(amrT)))

setdiff(dukeSamples$Sample, colnames(amrT))
setdiff(colnames(amrT), dukeSamples$Sample)

#RGI
rgiT<-read.delim("CountsTables/RGI_counts.tsv", sep = "\t", header = T, row.names = 1)
colnames(rgiT)<-gsub("pre", "PRE", colnames(rgiT), ignore.case = T) #switch all pre to cap
colnames(rgiT)<-sapply(str_split(colnames(rgiT), "_", n = 4 ), `[`, 1) #get rid of useless info
colnames(rgiT)<-gsub("\\.", "", colnames(rgiT))#get rid of random .s

length(intersect(dukeSamples$Sample, colnames(rgiT)))

setdiff(dukeSamples$Sample, colnames(rgiT))
setdiff(colnames(rgiT), dukeSamples$Sample)

#vsearch
vsearchT<-read.delim("CountsTables/vsearch_counts.tsv", sep = "\t", header = T, row.names = 1)
colnames(vsearchT)<-gsub("pre", "PRE", colnames(vsearchT), ignore.case = T) #switch all pre to cap
colnames(vsearchT)<-sapply(str_split(colnames(vsearchT), "_", n = 4 ), `[`, 1) #get rid of useless info
colnames(vsearchT)<-gsub("\\.", "", colnames(vsearchT))#get rid of random .s

length(intersect(dukeSamples$Sample, colnames(vsearchT)))

setdiff(dukeSamples$Sample, colnames(vsearchT))
setdiff(colnames(vsearchT), dukeSamples$Sample)
