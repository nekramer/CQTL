#!/usr/bin/R
library(tidyverse)
library(qvalue)
source("scripts/utils.R")

args <- commandArgs(trailingOnly = TRUE)
qtlData <- readQTLtools_perm(args[1])
geneInfo <- read_delim(args[2], delim = "\t", 
                       col_select = c("gene_id", "gene_name"))
outputFile <- args[3]

# Filter for non-NA and add column with q-values and column with FDR values
qtlData <- qtlData[!is.na(qtlData$adj_beta_val),] %>% 
  mutate("qval" = qvalue(adj_beta_pval)$qvalue) %>%
  mutate("FDR" = p.adjust(adj_beta_pval, method = "BH")) %>%
  # Join with geneInfo to get gene_name
  left_join(geneInfo) %>%
  # Write to file
  write_csv(file = outputFile)
