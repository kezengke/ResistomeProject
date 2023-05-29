#different antibiotic exposure windows metadata
rm(list = ls())

#exposure time frame
A_expo<-read.csv("AntibioticsExposureNewGroup.csv", check.names = F, header = T)
#metadata
metaData<-read.csv("metaWithBins.csv", header = T, row.names = 1, check.names = F)

IDtype<-unique(A_expo$Patient) #patient ID
ABtype<-unique(A_expo$Antibiotic) #antibiotic types

makeMeta <- function(days, metaData, ABtype, A_expo) {
  meta<-metaData
  for (i in 1:length(ABtype)) {
    AB<-vector()
    checkABframe<-A_expo[A_expo$Antibiotic == ABtype[i], ]
    for (j in 1:nrow(metaData)) {
      AB[j]<-sum(metaData[j, 2] > checkABframe[checkABframe$Patient == metaData[j, 1], 3] & 
                   metaData[j, 2] < (checkABframe[checkABframe$Patient == metaData[j, 1], 4]+days))
    }
    new_AB<-ifelse(AB > 0, "Yes", "No")
    cname<-ABtype[i]
    meta[, cname]<-new_AB
  }
  
  return(meta)
}

meta7days<-makeMeta(7, metaData, ABtype, A_expo)
meta10days<-makeMeta(10, metaData, ABtype, A_expo)
meta30days<-makeMeta(30, metaData, ABtype, A_expo)
meta60days<-makeMeta(60, metaData, ABtype, A_expo)
meta90days<-makeMeta(90, metaData, ABtype, A_expo)

write.csv(meta7days, "ABexposure/ABmeta7days.csv")
write.csv(meta10days, "ABexposure/ABmeta10days.csv")
write.csv(meta30days, "ABexposure/ABmeta30days.csv")
write.csv(meta60days, "ABexposure/ABmeta60days.csv")
write.csv(meta90days, "ABexposure/ABmeta90days.csv")


