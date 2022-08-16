#!/usr/bin/bash Rscript
library(data.table)
library(DESeq2)

args = commandArgs(trailingOnly = TRUE)
# args[1]: dds results from differential analysis
# args[2]: rsids

# Load in dds rda
load(args[1])

# Read in rsids
rsids <- fread(args[2], data.table = FALSE)$rsid

# Update rownames with rsids
rownames(allelic_imbalance_analysis) <- rsids

# Save allelic imbalance object
save(allelic_imbalance_analysis, file = "output/AI/differentialAllelicImbalance.rda")