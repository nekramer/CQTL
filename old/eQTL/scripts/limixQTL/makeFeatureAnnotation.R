#!/usr/bin/R
library(readr)
args <- commandArgs(trailingOnly = TRUE)
quantFile <- args[1]
outFile <- args[2]

# Read in quant file
CPMadjTMM_invNorm <- read_delim(quantFile, delim = "\t")


# Grab gene information
feature_annotation <- data.frame(
  "feature_id" = CPMadjTMM_invNorm$gene_id,
  "chromosome" = CPMadjTMM_invNorm$`#chr`,
  "start" = CPMadjTMM_invNorm$start,
  "end" = CPMadjTMM_invNorm$end,
  "ensembl_gene_id" = CPMadjTMM_invNorm$gene_id,
  "feature_strand" = CPMadjTMM_invNorm$strand,
  "gene_name" = CPMadjTMM_invNorm$gene_name
)

feature_annotation$chromosome <- gsub("chr", "",
                                      feature_annotation$chromosome)

# Write to file
write_delim(feature_annotation,
            file = outFile,
            delim = "\t")