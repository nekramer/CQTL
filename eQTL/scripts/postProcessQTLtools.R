#!/usr/bin/R
library(tidyverse)
source("scripts/utils.R")

get_eGene_variants <- function(eGene, permData, nomData){
  # Grab nominal p-value threshold for eGene
  nom_threshold <- permData %>% 
    filter(gene_id == eGene) %>% 
    pull(pval_nominal_threshold)
  
  # Filter nominal data for eGene and grab variants that satisfy nominal threshold
  eQTLs <- nomData %>% filter(gene_id == eGene) %>% filter(nom_pval <= nom_threshold)
  
  return(eQTLs)
  
}

args <- commandArgs(trailingOnly = TRUE)
threshold <- args[3]

geneInfo <- read_delim(args[4], delim = "\t", 
                       col_select = c("gene_id", "gene_name"))
outputFile <- args[5]

# Read in data from permutation pass, filter for non-NA values,
# add BH-adjusted p-value, and filter for significant eGenes based on threshold
permData <- readQTLtools_perm(args[2])
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

# Read in nominal data and grab variants that meet the pval_nominal_threshold
nomData <- readQTLtools_nom(args[1])
sig_nomData <- lapply(permData$gene_id, 
                      get_eGene_variants, 
                      permData = permData,
                      nomData = nomData) %>%
  bind_rows() %>%
  left_join(geneInfo) %>%
  write_csv(file = outputFile)
