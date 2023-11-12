#make pathway abundance counts table
rm(list = ls())
library(tidyverse)

file_pattern <- "*_pathabundance.tsv"
file_list <- list.files(path = "DukeHumann", pattern = file_pattern, recursive = TRUE, full.names = TRUE)

read_pathabundance <- function(file_name) {
  data <- read_tsv(file_name, col_types = cols(
    pathway = col_character(), 
    abundance = col_double()
  ),
  col_names = FALSE, skip = 1) %>%
    rename(pathway = X1, abundance = X2) %>%
    # Remove rows with UMAPPED/UNINTEGRATED
    filter(!grepl("UNMAPPED|UNINTEGRATED", pathway)) %>%
    mutate(sample = tools::file_path_sans_ext(basename(file_name)))
  
  return(data)
}

# Read each file and combine them into one data frame
pathabundance_data <- file_list %>% map_dfr(read_pathabundance)

# Pivot to wide format to create the counts table
counts_table <- pathabundance_data %>%
  pivot_wider(names_from = sample, values_from = abundance, values_fill = list(abundance = 0))


# Write the final counts table to a file
write_tsv(counts_table, "CountsTables/pathabundance_counts.tsv")
