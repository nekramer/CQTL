#!/usr/bin/R
library(tidyverse)
library(janitor)

args <- commandArgs(trailingOnly = TRUE)
prefix <- args[4]
Nk <- args[5]
validPEER <- read_csv(args[6], col_names = FALSE) %>% pull(X1)
batch <- as.logical(args[7])
numgenoPC <- as.numeric(args[9])

if (Nk %in% validPEER){
  # Read in PEER factors
  peer <- read_csv(args[1])
  
  # Read in donorInfo  
  donorInfo <- read_csv(args[2]) %>% dplyr::select(Donor, Sex)
  
  if (batch == TRUE){
    donorInfo <- read_csv(args[3]) %>% 
      distinct(Sample, .keep_all = TRUE) %>%
      filter(Condition == prefix) %>%
      arrange(match(Donor, peer$Donor)) %>%
      dplyr::select(Donor, RNAextractionKitBatch) %>%
      right_join(donorInfo)
  }
  
  # Join PEER factors with donorInfo to get Sex covariate
  covariates <- peer %>% left_join(donorInfo) %>% 
    # Transpose for QTLtools
    t() %>%
    as.data.frame() %>%
    # Rename rows and columns
    row_to_names(row_number = which(rownames(.) == "Donor"), 
                 remove_rows_above = FALSE) %>%
    rownames_to_column(var = "covariate")
  
  # Read in geno PCs
  genoPCs <- read_delim(args[8], delim = " ") %>% 
    dplyr::select(-SampleID) %>%
    # Put donor columns in same order as RNA PCs
    relocate(all_of(peer$Donor)) %>%
    mutate("covariate" = paste0("geno_PC", 1:nrow(.))) %>%
    relocate(covariate) %>%
    # Keep different number of PC's depending on condition
    filter(covariate %in% paste0("geno_PC", 1:numgenoPC))
  
  # Join all covariates
  covariates <- rbind(covariates, genoPCs)
  
  # Write to tab-delimited file
  if (batch == TRUE){
    write_delim(covariates, file = paste0("output/covar/", prefix, "_PEER_k", Nk,"_genoPC_batch_covar.txt"), delim = "\t")
  } else {
    write_delim(covariates, file = paste0("output/covar/", prefix, "_PEER_k", Nk,"_genoPC_covar.txt"), delim = "\t")
  }
  
  
} else {
  
  if (batch == TRUE){
    file.create(paste0("output/covar/", prefix, "_PEER_k", Nk,"_genoPC_batch_covar.txt"))
  } else {
    file.create(paste0("output/covar/", prefix, "_PEER_k", Nk,"_genoPC_covar.txt"))
  }
  
}