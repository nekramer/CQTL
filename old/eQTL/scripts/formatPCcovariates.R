#!/usr/bin/R
library(tidyverse)
library(janitor)

args <- commandArgs(trailingOnly = TRUE)
numPC <- as.numeric(args[4])
prefix <- args[5]
batch <- as.logical(args[6])

# Read in PCs 
PCs <- read_csv(args[1]) %>% 
  dplyr::rename("Donor" = 1) %>%
  # Keep different number of PC's depending on condition
  dplyr::select(Donor, num_range("PC", 1:numPC))
  
  
donorInfo <- read_csv(args[2]) %>% dplyr::select(Donor, Sex)

if (batch == TRUE){
  donorInfo <- read_csv(args[3]) %>% 
    distinct(Sample, .keep_all = TRUE) %>%
    filter(Condition == prefix) %>%
    arrange(match(Donor, PCs$Donor)) %>%
    dplyr::select(Donor, RNAextractionKitBatch) %>%
    right_join(donorInfo)
}

# Join PC's with donorInfo to get Sex covariate
covariates <- PCs %>% left_join(donorInfo) %>% 
  # Transpose for QTLtools
  t() %>%
  as.data.frame() %>%
  # Rename rows and columns
  rownames_to_column(var = "covariate") %>%
  row_to_names(row_number = 1)
colnames(covariates)[1] <- "covariate"

# Write to tab-delimited file
if (batch == TRUE){
  write_delim(covariates, file = paste0("output/covar/", prefix, "_PC_batch_covar.txt"),
              delim = "\t")
} else {
  write_delim(covariates, file = paste0("output/covar/", prefix, "_PC_covar.txt"),
              delim = "\t")
}

  


