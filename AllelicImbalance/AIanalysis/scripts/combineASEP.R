library(dplyr)
library(readr)

# Load in the the asep objects, add chrom column, and add to list
asep_list <- list()
for (chr in 1:22){
  load(paste0("data/chr", chr, "_ASEP.rda"))
  chrASEP_name <- paste0("chr", chr, "_asep")
  asep_results <- as.data.frame(asep_results)
  asep_results$chrom <- paste0("chr", chr)
  
  # Add to asep_list with name
  asep_list[[chrASEP_name]] <- asep_results
}

# Bind all
asep_all <- bind_rows(asep_list)

# Write to file
write_csv(asep_all, file = "data/ASEP_results.csv")