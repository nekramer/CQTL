library(dplyr)
library(tibble)

# args[1] is the chromosome
args = commandArgs(trailingOnly = TRUE)

# Read in lifted over GRCh38 file, sorting one by the physical position and
# the other by the genomic position, assigning each rowid as a column value
gm38_physical <- read.table(paste0("genetic_map_GRCh38_chr", args[1], ".txt"), 
                   header = TRUE) %>% arrange(PhysicalPosition.bp.) %>% rowid_to_column()
gm38_genomic <- read.table(paste0("genetic_map_GRCh38_chr", args[1], ".txt"), 
                            header = TRUE) %>% arrange(GenomicPosition.cM.) %>% rowid_to_column()

# Semi-join to get all rows that are the same order in same index position
gm38_sorted <- semi_join(gm38_physical, gm38_genomic)
# Remove rowid column
gm38_sorted <- gm38_sorted[,-1]
 
# Write new file
write.table(gm38_sorted, file = paste0("genetic_map_GRCh38_chr", args[1], "_sorted.txt"),
            quote = FALSE, row.names = FALSE, col.names = TRUE)
