#!/usr/bin/env Rscript

library(data.table)
library(dplyr)

args = commandArgs(trailingOnly = TRUE)

# Read in data bim file
data <- fread(args[1], data.table = FALSE)

# Read in ref bim file
ref <- fread(args[2], data.table = FALSE)

# Join data and ref by variant location
joined <- inner_join(data, ref, by = c("V1", "V4"))

# Compare alleles (V5.x and V5.y & V6.x and V6.y)
missnps <- joined[which(joined$V5.x != joined$V5.y | joined$V6.x != joined$V6.y),]

# Write lists of snps to txt files to remove from data/reference
write.table(missnps$V2.x, file = "output/ancestry/data_missnp.txt", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(missnps$V2.y, file = "output/ancestry/ref_missnp.txt", quote = FALSE, col.names = FALSE, row.names = FALSE)