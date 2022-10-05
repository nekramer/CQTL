#!/usr/bin/R
library(tidyverse)
library(janitor)

args <- commandArgs(trailingOnly = TRUE)
prefix <- args[3]

# Read in PCs and split Sample column into Donor column
PCs <- read_csv(args[1]) %>% 
  dplyr::rename("Sample" = 1) %>% 
  separate(Sample, into = c(NA, "Donor", NA, NA, NA, NA), sep = "_", remove = FALSE)
donorInfo <- read_csv(args[2]) %>% dplyr::select(Donor, Sex)

# Join PC's with donorInfo to get Sex covariate
covariates <- PCs %>% left_join(donorInfo) %>% 
  dplyr::select(-Sample) %>% 
  # Transpose for QTLtools
  t() %>%
  as.data.frame() %>%
  # Rename rows and columns
  rownames_to_column(var = "covariate") %>%
  row_to_names(row_number = 1)
colnames(covariates)[1] <- "covariate"

# Write to tab-delimited file
write_delim(covariates, file = paste0("output/covar/", prefix, "_PCcovar.txt"),
            delim = "\t")
  


