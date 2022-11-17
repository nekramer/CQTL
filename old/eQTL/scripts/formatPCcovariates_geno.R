#!/usr/bin/R
library(tidyverse)
library(janitor)

args <- commandArgs(trailingOnly = TRUE)


numPC <- as.numeric(args[4])
numgenoPC <- as.numeric(args[6])
batch <- as.logical(args[7])
prefix <- args[8]


# Read in RNA PCs 
PCs <- read_csv(args[1]) %>% 
  dplyr::rename("Donor" = 1) %>%
  # Keep different number of PC's depending on condition
  dplyr::select(Donor, num_range("PC", 1:numPC))

# Read in geno PCs
genoPCs <- read_delim(args[5], delim = " ") %>% 
  dplyr::select(-SampleID) %>%
  # Put donor columns in same order as RNA PCs
  relocate(all_of(PCs$Donor)) %>%
  mutate("covariate" = paste0("geno_PC", 1:nrow(.))) %>%
  relocate(covariate) %>%
  # Keep different number of PC's depending on condition
  filter(covariate %in% paste0("geno_PC", 1:numgenoPC))

donorInfo <- read_csv(args[2]) %>% dplyr::select(Donor, Sex)

if (batch == TRUE){
  donorInfo <- read_csv(args[3]) %>% 
    distinct(Sample, .keep_all = TRUE) %>%
    filter(Condition == prefix) %>%
    arrange(match(Donor, PCs$Donor)) %>%
    dplyr::select(Donor, RNAextractionKitBatch) %>%
    right_join(donorInfo)
}

# Join all covariates
covariates <- PCs %>% left_join(donorInfo) %>% 
  # Transpose
  t() %>%
  as.data.frame() %>%
  # Rename rows and columns
  rownames_to_column(var = "covariate") %>%
  row_to_names(row_number = 1)
colnames(covariates)[1] <- "covariate"

covariates <- rbind(covariates, genoPCs)

# Write to tab-delimited file
if (batch == TRUE){
  write_delim(covariates, file = paste0("output/covar/", prefix, 
                                        "_PC_genoPC_batch_covar.txt"),
              delim = "\t")
} else {
  write_delim(covariates, file = paste0("output/covar/", prefix, 
                                        "_PC_genoPC_covar.txt"),
              delim = "\t")
}




