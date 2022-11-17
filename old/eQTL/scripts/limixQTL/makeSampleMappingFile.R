#!/usr/bin/R
library(readr)
library(dplyr)
args <- commandArgs(trailingOnly = TRUE)
outFile <- args[2]

donorInfo <- read_csv(args[1]) %>% pull(Donor)

sampleMappingFile <- data.frame(
  "geno_individual_id" = donorInfo,
  "phenotype_sample_id" = donorInfo
)

write_delim(sampleMappingFile, file = outFile, delim = "\t")