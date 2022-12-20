#!/usr/bin/R
library(tidyverse)
library(janitor)

args <- commandArgs(trailingOnly = TRUE)

# Read in all PEER factors, rename Donor to Sample and extract Donor
peer <- read_csv(args[1]) %>%
  dplyr::rename(Sample = Donor) %>%
  separate(Sample, into = c(NA, "Donor", NA, "Condition", NA, NA), sep = "_",
           remove = FALSE)
  
# Read in donorInfo to get sex (age?)
donorInfo <- read_csv(args[2]) %>% dplyr::select(Donor, Sex)


# Read in DNA samplesheet
dnaBatches <- c("Batch", "DNAReagentBatch")[c(as.logical(args[4]), 
                                              as.logical(args[5]))]
dnaInfo <- read_csv(args[3]) %>% dplyr::select(Donor, all_of(dnaBatches))

# Read in RNA samplesheet
rnaBatches <- c("RNAextractionKitBatch", "SequencingBatch")[c(as.logical(args[7]),
                                                              as.logical(args[8]))]
rnaInfo <- read_csv(args[6]) %>% 
  filter(Sample %in% peer$Sample) %>%
  distinct(Sample, .keep_all = TRUE) %>%
  dplyr::select(Sample, Donor, Condition, all_of(rnaBatches))

# Join covariates thus far
covariates <- peer %>%
  full_join(donorInfo) %>%
  full_join(dnaInfo) %>%
  full_join(rnaInfo)

# Read in genoPCs and reformat to join with covariates
numgenoPC <- as.numeric(args[10])
genoPCs <- read_delim(args[9], delim = " ") %>%
  dplyr::select(-SampleID) %>%
  mutate("covariate" = paste0("geno_PC", 1:nrow(.))) %>%
  # Filter for numgenoPC
  filter(covariate %in% paste0("geno_PC", 1:numgenoPC)) %>%
  t() %>%
  as.data.frame() %>%
  # Rename rows and columns
  row_to_names(row_number = which(rownames(.) == "covariate"), 
               remove_rows_above = FALSE) %>%
  rownames_to_column(var = "Donor")


# Join covariates with genoPCs
covariates <- covariates %>% 
  full_join(genoPCs) %>%
  # Remove Condition and Donor column now that everything is joined
  dplyr::select(-Donor, -Condition) %>%
  relocate(Sample)


# Write to file
write_csv(covariates, file = args[11])