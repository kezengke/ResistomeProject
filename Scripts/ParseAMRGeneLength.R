rm(list = ls())

startEnd<-read.table("ReferenceGeneCatalog.txt", sep = "\t", header = T, fill = T, quote = "")
startEnd$gene_length<-startEnd$refseq_stop - startEnd$refseq_start

write.csv(startEnd, "AMRGeneLength.csv")
