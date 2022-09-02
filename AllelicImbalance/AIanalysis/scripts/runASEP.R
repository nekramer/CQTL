library(ASEP)
library(readr)
library(dplyr)
library(parallel)

args = commandArgs(trailingOnly = TRUE)
# args[1] is the phased file
# args[2] is the chromosome

data_phased <- read_csv(args[1]) %>%
  # Remove SNPs that don't have a gene
  filter(!is.na(gene)) %>%
  # Set "ref_condition"
  mutate(ref_condition = "CTL") %>%
  # Group by gene
  group_by(gene) %>%
  # Filter out genes that don't have snps in both CTL and FNF
  filter(length(unique(group)) == 2) %>%
  ungroup() %>%
  # Now group by id
  group_by(gene, id) %>%
  # Filter out donors in each gene that don't have a CTL and FNF
  filter(length(unique(group)) == 2)

# Run ASEP
asep_results <- differential_ASE_detection(data_phased, 
                           phased = TRUE, 
                           varList = NULL,
                           parallel = TRUE,
                           n_core = 8)

# Write file
save(asep_results, file = paste0("data/chr", args[2], "_ASEP.rda"))

#write_csv(asep_results, file = paste0("data/chr", args[2], "_ASEP.txt"))

