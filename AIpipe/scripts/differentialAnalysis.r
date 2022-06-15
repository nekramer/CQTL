#!/usr/bin/bash Rscript
library(data.table)
library(DESeq2)
library(stringr)

args = commandArgs(trailingOnly = TRUE)
# args[1]: alleleCounts matrix
# args[2]: colData for alleleCounts matrix
# args[3]: weight matrix

# Read in filtered alleleCounts matrix
alleleCountsMatrix <- read.table(args[1])

# Convert rownames to separate dataframe of chrom and position
positions <- data.frame(do.call(rbind, str_split(rownames(alleleCountsMatrix), ":")))
positions <- positions[,c(1,2)]
colnames(positions) <- c("chr", "pos")
write.table(positions, file = "output/AI/AI_variant_positions.txt", quote = FALSE, row.names = FALSE, col.names = TRUE)

# Read in sampleData 
colData <- fread(args[2], data.table = FALSE)

# Read in weight matrix and set row names
weightMatrix <- fread(args[3], data.table = FALSE)
rownames(weightMatrix) <- weightMatrix$V1
weightMatrix <- as.matrix(weightMatrix[,-1])

# Initialize design for DESeq
design <- ~0 + Treatment:Donor + Treatment:Allele

# Add datasets, weights, and size factors
dds <- DESeqDataSetFromMatrix(alleleCountsMatrix, colData, design)
assays(dds)[["weights"]] <- weightMatrix
# No size factor normalization, set all to 1
sizeFactors(dds) <- rep(1, ncol(alleleCountsMatrix))
# Relevel with reference allele set first
dds$Allele <- relevel(dds$Allele, ref = "ref")

allelic_imbalance_analysis <- DESeq(dds)

# Save allelic imbalance object
save(allelic_imbalance_analysis, file = "output/AI/differentialAllelicImbalance_v1.rda")