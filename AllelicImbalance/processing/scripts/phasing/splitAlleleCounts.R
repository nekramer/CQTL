library(data.table)
library(tidyr)

writeSplitFile <- function(chromName){
  fwrite(chrom_alleleCounts_hets[[chromName]],
         file = paste0('../../output/AI/', chromName,'_alleleCounts_5hets.csv'),
         row.names = FALSE, quote = FALSE)
}


# Read in data
alleleCounts_hets <- tibble(fread('../../output/AI/alleleCounts_5hets.csv')) %>%
  # Split variantID column to get chrom, position, refAllele, and altAllele columns
  separate(variantID, c("chrom", "pos", "refAllele", "altAllele"))

# Split into separate chromosome dataframes
chrom_alleleCounts_hets <- split(alleleCounts_hets, alleleCounts_hets$chrom)

# Write to files
invisible(lapply(names(chrom_alleleCounts_hets), writeSplitFile))