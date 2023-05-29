rm(list = ls())
ABtypes<-read.csv("ABtypes.csv")
ABtypes$Antibiotic<-tolower(ABtypes$Antibiotic)
ABtypes$Antibiotic<-trimws(ABtypes$Antibiotic, which = "right")
sort(unique(ABtypes$Antibiotic))

write.csv(sort(unique(ABtypes$Antibiotic)), "AllABtypes.csv")
