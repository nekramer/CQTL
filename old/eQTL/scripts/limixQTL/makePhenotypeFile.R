#!/usr/bin/R
library(readr)
library(dplyr)
args <- commandArgs(trailingOnly = TRUE)
quantFile <- args[1]
outFile <- args[2]

# Read in quant file
CPMadjTMM_invNorm <- read_delim(quantFile, delim = "\t")

# Grab geneid and sample gene quant info
phenotype_file <- CPMadjTMM_invNorm %>% 
  dplyr::select(-`#chr`, -start, -end, -gene_name, -strand) %>%
  dplyr::rename("feature_id" = gene_id)

# Write to file
write_delim(phenotype_file,
            file = outFile,
            delim = "\t")