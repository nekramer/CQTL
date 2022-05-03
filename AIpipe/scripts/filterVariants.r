#!/usr/bin/env Rscript

library(data.table)
args = commandArgs(trailingOnly = TRUE)
# args[1]: alleleCounts matrix
# args[2]: variants to remove
# args[3]: RNAhets

# Read in alleleCounts matrix
alleleCounts <- read.table(args[1])
# Read in RNAhets
RNAhets <- fread(args[3], data.table = FALSE)

# Read in list of variants to remove
remove <- fread(args[2], data.table = FALSE)

# subset alleleCountsMatrix for variants that aren't in remove
alleleCounts_subset <- alleleCounts[!rownames(alleleCounts) %in% remove$variant,]
# subset RNAhets for variants that aren't in remove
RNAhets_subset <- RNAhets[!RNAhets$variant %in% remove$variant,]

# Write to files
write.table(alleleCounts_subset, "output/AI/alleleCountsMatrix_filtered.txt", quote = FALSE, row.names = TRUE, col.names = TRUE)
write.table(RNAhets_subset, "output/AI/RNAhets_filtered.txt", quote = FALSE, row.names = FALSE, col.names = TRUE)