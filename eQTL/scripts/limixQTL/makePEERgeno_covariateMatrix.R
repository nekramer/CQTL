#!/usr/bin/R
library(tidyverse)
library(janitor)
args <- commandArgs(trailingOnly = TRUE)
condition <- args[4]
Nk <- args[5]
numgenoPC <- as.numeric(args[7])

# Read in PEER factors
peer <- read_csv(args[1])

# Read in donorInfo  
donorInfo <- read_csv(args[2]) %>% dplyr::select(Donor, Sex)

# Get RNA batch information from samplesheet
donorInfo <- read_csv(args[3]) %>% 
  distinct(Sample, .keep_all = TRUE) %>%
  filter(Condition == condition) %>%
  arrange(match(Donor, peer$Donor)) %>%
  dplyr::select(Donor, RNAextractionKitBatch) %>%
  right_join(donorInfo)

# Join PEER factors with donorInfo
covariates <- peer %>% left_join(donorInfo)

# Read in geno PCs
genoPCs <- read_delim(args[6], delim = " ") %>% 
  dplyr::select(-SampleID) %>%
  # Put donor columns in same order as RNA PCs
  relocate(all_of(peer$Donor)) %>%
  mutate("covariate" = paste0("geno_PC", 1:nrow(.))) %>%
  relocate(covariate) %>%
  # Keep number of genoPCs specified by kneedle
  filter(covariate %in% paste0("geno_PC", 1:numgenoPC)) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "covariate") %>%
  row_to_names(row_number = 1) %>%
  dplyr::rename("Donor" = covariate)

# Join covariates with genoPCs
covariates <- covariates %>% left_join(genoPCs) %>%
  relocate(Donor) %>%
  dplyr::rename("sample_id" = Donor)

write_delim(covariates, file = paste0("output/limix_input/", condition, "_PEER_k", Nk,"_genoPC_batch_covariateMatrix.tsv"), delim = "\t")