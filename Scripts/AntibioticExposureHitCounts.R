rm(list = ls())

bracken7<-read.csv("ABexposure/bracken7.csv", header = T, row.names = 1)
bracken10<-read.csv("ABexposure/bracken10.csv", header = T, row.names = 1)
bracken30<-read.csv("ABexposure/bracken30.csv", header = T, row.names = 1)
bracken60<-read.csv("ABexposure/bracken60.csv", header = T, row.names = 1)
bracken90<-read.csv("ABexposure/bracken90.csv", header = T, row.names = 1)

amr7<-read.csv("ABexposure/amr7.csv", header = T, row.names = 1)
amr10<-read.csv("ABexposure/amr10.csv", header = T, row.names = 1)
amr30<-read.csv("ABexposure/amr30.csv", header = T, row.names = 1)
amr60<-read.csv("ABexposure/amr60.csv", header = T, row.names = 1)
amr90<-read.csv("ABexposure/amr90.csv", header = T, row.names = 1)

rgi7<-read.csv("ABexposure/rgi7.csv", header = T, row.names = 1)
rgi10<-read.csv("ABexposure/rgi10.csv", header = T, row.names = 1)
rgi30<-read.csv("ABexposure/rgi30.csv", header = T, row.names = 1)
rgi60<-read.csv("ABexposure/rgi60.csv", header = T, row.names = 1)
rgi90<-read.csv("ABexposure/rgi90.csv", header = T, row.names = 1)

vsearch7<-read.csv("ABexposure/vsearch7.csv", header = T, row.names = 1)
vsearch10<-read.csv("ABexposure/vsearch10.csv", header = T, row.names = 1)
vsearch30<-read.csv("ABexposure/vsearch30.csv", header = T, row.names = 1)
vsearch60<-read.csv("ABexposure/vsearch60.csv", header = T, row.names = 1)
vsearch90<-read.csv("ABexposure/vsearch90.csv", header = T, row.names = 1)

HitsT<-rbind(colSums(bracken7<0.05), colSums(bracken10<0.05), colSums(bracken30<0.05), colSums(bracken60<0.05), colSums(bracken90<0.05),
             colSums(amr7<0.05), colSums(amr10<0.05), colSums(amr30<0.05), colSums(amr60<0.05), colSums(amr90<0.05),
             colSums(rgi7<0.05), colSums(rgi10<0.05), colSums(rgi30<0.05), colSums(rgi60<0.05), colSums(rgi90<0.05),
             colSums(vsearch7<0.05), colSums(vsearch10<0.05), colSums(vsearch30<0.05), colSums(vsearch60<0.05), colSums(vsearch90<0.05))
colnames(HitsT)<-colnames(bracken7)
rownames(HitsT)<-c("brackenD7", "brackenD10", "brackenD30", "brackenD60", "brackenD90",
                   "amrD7", "amrD10", "amrD30", "amrD60", "amrD90",
                   "rgiD7", "rgiD10", "rgiD30", "rgiD60", "rgiD90",
                   "vsearchD7", "vsearchD10", "vsearchD30", "vsearchD60", "vsearchD90")

write.csv(HitsT, "AntibioticExposureHitCounts.csv")

#adjusted
rm(list = ls())

bracken7<-read.csv("ABexposure/bracken7adj.csv", header = T, row.names = 1)
bracken10<-read.csv("ABexposure/bracken10adj.csv", header = T, row.names = 1)
bracken30<-read.csv("ABexposure/bracken30adj.csv", header = T, row.names = 1)
bracken60<-read.csv("ABexposure/bracken60adj.csv", header = T, row.names = 1)
bracken90<-read.csv("ABexposure/bracken90adj.csv", header = T, row.names = 1)

amr7<-read.csv("ABexposure/amr7adj.csv", header = T, row.names = 1)
amr10<-read.csv("ABexposure/amr10adj.csv", header = T, row.names = 1)
amr30<-read.csv("ABexposure/amr30adj.csv", header = T, row.names = 1)
amr60<-read.csv("ABexposure/amr60adj.csv", header = T, row.names = 1)
amr90<-read.csv("ABexposure/amr90adj.csv", header = T, row.names = 1)

rgi7<-read.csv("ABexposure/rgi7adj.csv", header = T, row.names = 1)
rgi10<-read.csv("ABexposure/rgi10adj.csv", header = T, row.names = 1)
rgi30<-read.csv("ABexposure/rgi30adj.csv", header = T, row.names = 1)
rgi60<-read.csv("ABexposure/rgi60adj.csv", header = T, row.names = 1)
rgi90<-read.csv("ABexposure/rgi90adj.csv", header = T, row.names = 1)

vsearch7<-read.csv("ABexposure/vsearch7adj.csv", header = T, row.names = 1)
vsearch10<-read.csv("ABexposure/vsearch10adj.csv", header = T, row.names = 1)
vsearch30<-read.csv("ABexposure/vsearch30adj.csv", header = T, row.names = 1)
vsearch60<-read.csv("ABexposure/vsearch60adj.csv", header = T, row.names = 1)
vsearch90<-read.csv("ABexposure/vsearch90adj.csv", header = T, row.names = 1)

HitsT<-rbind(colSums(bracken7<0.05), colSums(bracken10<0.05), colSums(bracken30<0.05), colSums(bracken60<0.05), colSums(bracken90<0.05),
             colSums(amr7<0.05), colSums(amr10<0.05), colSums(amr30<0.05), colSums(amr60<0.05), colSums(amr90<0.05),
             colSums(rgi7<0.05), colSums(rgi10<0.05), colSums(rgi30<0.05), colSums(rgi60<0.05), colSums(rgi90<0.05),
             colSums(vsearch7<0.05), colSums(vsearch10<0.05), colSums(vsearch30<0.05), colSums(vsearch60<0.05), colSums(vsearch90<0.05))
colnames(HitsT)<-colnames(bracken7)
rownames(HitsT)<-c("brackenD7adj", "brackenD10adj", "brackenD30adj", "brackenD60adj", "brackenD90adj",
                   "amrD7adj", "amrD10adj", "amrD30adj", "amrD60adj", "amrD90adj",
                   "rgiD7adj", "rgiD10adj", "rgiD30adj", "rgiD60adj", "rgiD90adj",
                   "vsearchD7adj", "vsearchD10adj", "vsearchD30adj", "vsearchD60adj", "vsearchD90adj")

write.csv(HitsT, "AntibioticExposureHitCountsAdj.csv")
