#!/usr/bin/R
library(tidyverse)
library(janitor)

args <- commandArgs(trailingOnly = TRUE)

# Read in PEER factors and subset for numPEER
if (file.exists(args[2])){
  numPEER <- read_csv(args[2], col_names = FALSE) %>% pull(X1)
  peer <- read_csv(args[1]) %>%
    dplyr::select(paste0("PEER", 1:numPEER), Donor)
} else {
  peer <- read_csv(args[1])
}

# Read in donorInfo to get sex (age?)
donorInfo <- read_csv(args[3]) %>% dplyr::select(Donor, Sex)


# Read in DNA samplesheet
dnaBatches <- c("Batch", "DNAReagentBatch")[c(as.logical(args[5]), 
                                              as.logical(args[6]))]
dnaInfo <- read_csv(args[4]) %>% dplyr::select(Donor, all_of(dnaBatches))

# Read in RNA samplesheet
rnaBatches <- c("RNAextractionKitBatch", "SequencingBatch")[c(as.logical(args[8]),
                                                              as.logical(args[9]))]
condition <- args[10]
rnaInfo <- read_csv(args[7]) %>% 
  distinct(Sample, .keep_all = TRUE) %>%
  filter(Condition == condition) %>%
  dplyr::select(Donor, all_of(rnaBatches))

# Join covariates thus far
covariates <- peer %>%
  full_join(donorInfo) %>%
  full_join(dnaInfo) %>%
  full_join(rnaInfo) %>%
  # Transpose
  t() %>%
  as.data.frame() %>%
  # Rename rows and columns
  row_to_names(row_number = which(rownames(.) == "Donor"), 
               remove_rows_above = FALSE) %>%
  rownames_to_column(var = "covariate")

# Read in genoPCs
numgenoPC <- as.numeric(args[12])
genoPCs <- read_delim(args[11], delim = " ") %>%
  dplyr::select(-SampleID) %>%
  mutate("covariate" = paste0("geno_PC", 1:nrow(.))) %>%
  # Put donor columns in same order of covariates
  relocate(all_of(colnames(covariates))) %>%
  # Filter for numgenoPC
  filter(covariate %in% paste0("geno_PC", 1:numgenoPC))

# Join covariates with genoPCs
covariates <- rbind(covariates, genoPCs)

# Write to tab-delimited file
write_delim(covariates, file = args[13], delim = "\t")
  