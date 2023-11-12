#normalizing counts tables
rm(list = ls())
library(stringr)

#counts tables
brackenT<-read.csv("CountsTables/brackenProcessed.csv", header = T, row.names = 1, check.names = F)
genusT<-read.csv("CountsTables/genusProcessed.csv", header = T, row.names = 1, check.names = F)
phylumT<-read.csv("CountsTables/phylumProcessed.csv", header = T, row.names = 1, check.names = F)

amrT<-read.csv("CountsTables/amrProcessed.csv", header = T, row.names = 1, check.names = F)
rgiT<-read.csv("CountsTables/rgiProcessed.csv", header = T, row.names = 1, check.names = F)
vsearchT<-read.csv("CountsTables/vsearchProcessed.csv", header = T, row.names = 1, check.names = F)

pathT<-read.csv("CountsTables/pathwayProcessed.csv", header = T, row.names = 1, check.names = F)

amrGeneLength<-read.csv("AMRGeneLength.csv", header = T)
matched_indices <- sapply(rownames(amrT), function(name) {
  which(apply(amrGeneLength, 1, function(row) name %in% row))[1]
})
matched_indices <- matched_indices[!is.na(matched_indices)]
amrGeneLength<-amrGeneLength[matched_indices, ]

rgiGeneLength<-read.csv("CARD_geneLength.csv", header = T, row.names = 1)
rgiGeneLength<-rgiGeneLength[match(rownames(rgiT), rgiGeneLength$ARO_name), , drop = F]

vsearchGeneLength<-read.csv("VsearchGeneLength.csv", header = T)
vsearchGeneLength$Header<-sapply(str_split(sapply(str_split(vsearchGeneLength$Header, "\\|", n = 6), `[`, 6), " \\[", n = 2), `[`, 1)
vsearchGeneLength<-vsearchGeneLength[match(rownames(vsearchT), vsearchGeneLength$Header), , drop = F]

microbNorm <- function(counts) {
  #Normalization
  n<-colSums(counts)
  sumx<-sum(counts)
  for (i in 1:ncol(counts)) {
    counts[,i]<-counts[,i]/n[i]
  }
  counts<-log10(counts*(sumx/ncol(counts))+1)
  
  counts<-data.frame(counts, check.names = F)
  return(counts)
}

geneNorm <- function(counts, gene_lengths) {
  # Normalize by gene length to get RPK
  RPK <- t(t(counts) / (gene_lengths / 1000))
  
  # Get the total number of reads (in millions)
  total_reads <- colSums(counts) / 10^6
  
  # Normalize RPK by total reads to get RPKM
  RPKM <- t(t(RPK) / total_reads)
  
  RPKM <- data.frame(RPKM, check.names = F)
  return(RPKM)
}

speciesNorm<-microbNorm(brackenT)
genusNorm<-microbNorm(genusT)
phylumNorm<-microbNorm(phylumT)

amrNorm<-geneNorm(amrT, amrGeneLength$gene_length)
rgiNorm<-geneNorm(rgiT, rgiGeneLength$gene_length)
vsearchNorm<-geneNorm(vsearchT, vsearchGeneLength$Length)

pathNorm<-microbNorm(pathT)

write.csv(speciesNorm, "CountsTables/brackenNormalized.csv")
write.csv(genusNorm, "CountsTables/genusNormalized.csv")
write.csv(phylumNorm, "CountsTables/phylumNormalized.csv")
write.csv(amrNorm, "CountsTables/amrNormalized.csv")
write.csv(rgiNorm, "CountsTables/rgiNormalized.csv")
write.csv(vsearchNorm, "CountsTables/vsearchNormalized.csv")
write.csv(pathNorm, "CountsTables/pathNormalized.csv")

