#!/usr/bin/R
library(tidyverse)
source("scripts/utils.R")

args <- commandArgs(trailingOnly = TRUE)

threshold <- as.numeric(args[2])
outFile <- args[3]
# Read in data from permutation pass, filter for non-NA values,
# add BH-adjusted p-value, and filter for significant eGenes based on threshold
permData <- readQTLtools_perm(args[1])

permData <- permData[!is.na(permData$adj_beta_pval),] %>%
  mutate("adj_FDR" = p.adjust(adj_beta_pval, method = "BH"))

# Smallest p-value above FDR threshold
upper_bound <- permData %>% filter(adj_FDR > threshold) %>% arrange(adj_beta_pval) %>%
  dplyr::slice(1) %>% pull(adj_beta_pval)

# Largest p-value below FDR threshold
lower_bound <- permData %>% filter(adj_FDR <= threshold) %>% 
  arrange(desc(adj_beta_pval)) %>% dplyr::slice(1) %>%
  pull(adj_beta_pval)

# Calculate nominal p-value threshold for each gene
pthreshold <- (lower_bound + upper_bound)/2
permData <- permData %>% mutate("pval_nominal_threshold" =
                                  qbeta(pthreshold, bml1, bml2, ncp = 0, lower.tail = TRUE, log.p = FALSE))

# Write gene id and pval_nominal_threshold to outfile
permData %>% dplyr::select(gene_id, pval_nominal_threshold) %>%
  write_csv(file = outFile)

