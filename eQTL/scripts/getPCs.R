#!/usr/bin/R
library(tidyverse)
library(janitor)

args <- commandArgs(trailingOnly = TRUE)
samplesheet <- args[1]
normQuant <- args[2]
prefix <- args[3]

samplesheet <- read_csv(samplesheet) %>%
  # Get distinct samples
  distinct(Sample, .keep_all = TRUE)

CPMadjTMM_invNorm <- read_delim(normQuant, delim = "\t")

CPMadjTMM_invNorm_info <- CPMadjTMM_invNorm %>% 
  dplyr::select(-`#chr`, -start, -end, -strand, -gene_name) %>%
  t() %>% as.data.frame() %>% row_to_names(row_number = 1) %>%
  rownames_to_column("Sample") %>%
  left_join(samplesheet) 

# Calculate variances to remove columns where variance is 0 (can't do PCA)
variances <- CPMadjTMM_invNorm %>% 
  dplyr::select(-gene_id, -`#chr`, -start, -end, -strand, -gene_name) %>% 
  t() %>% apply(2, var, na.rm = TRUE)

if (length(which(variances == 0) > 0)){
  
  CPMadjTMM_invNorm_pca <- CPMadjTMM_invNorm %>% 
    dplyr::select(-gene_id, -`#chr`, -start, -end, -strand, -gene_name) %>% 
    t() %>% as.data.frame() %>% dplyr::select(-which(variances == 0)) %>%
    prcomp(center = TRUE, scale = TRUE)
  
} else {
  
  CPMadjTMM_invNorm_pca <- CPMadjTMM_invNorm %>% 
    dplyr::select(-gene_id, -`#chr`, -start, -end, -strand, -gene_name) %>% 
    t() %>%
    prcomp(center = TRUE, scale = TRUE)
}

# Bind first 10 PC's to rest of info
CPMadjTMM_invNorm_info <- bind_cols(CPMadjTMM_invNorm_info, CPMadjTMM_invNorm_pca$x[,1:10])

# Write file
write_csv(bind_cols(CPMadjTMM_invNorm_info$Sample, CPMadjTMM_invNorm_pca$x[,1:10]), file = paste0("output/covar/", prefix, "_PC.csv"))
