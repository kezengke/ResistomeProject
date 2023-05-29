rm(list = ls())
newGroup<-read.csv("AllABtypes_group.csv")
ABdays<-read.csv("AntibioticsExposure.csv")
ABnames<-read.csv("ABtypes.csv")

AB<-data.frame(ABnames$PATIENT, ABnames$Antibiotic, ABdays$DaysStart, ABdays$DaysEnd)
colnames(AB)<-c("Patient", "Antibiotic", "DaysStart", "DaysEnd")

AB$Antibiotic<-tolower(AB$Antibiotic)
AB$Antibiotic<-trimws(AB$Antibiotic, which = "right")

AB$Antibiotic<-newGroup$Antibiotic[match(AB$Antibiotic, newGroup$AllABtypes)]

write.csv(AB, "AntibioticsExposureNewGroup.csv", row.names = F)
