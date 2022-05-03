#!/usr/bin/env Rscript

library(data.table)
library(stringr)
args = commandArgs(trailingOnly = TRUE)
# args[1]: donors in RNA samplesheet
# args[2]: file that has corresponding RNA/genotype donor names

donors <- unlist(str_split(args[1], ","))
donorNames <- fread(args[2], data.table = FALSE)

# Subset donorNames for which donors are found in the samplesheet
donorNames <- donorNames[which(donorNames$rna %in% donors),]

# Write remaining 'vcf' column to file
write.table(donorNames$vcf, "output/vcf/AIsamples.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)