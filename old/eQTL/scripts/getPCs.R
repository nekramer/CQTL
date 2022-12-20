#!/usr/bin/R
library(tidyverse)
library(janitor)

args <- commandArgs(trailingOnly = TRUE)
samplesheet <- args[1]
normQuant <- args[2]
prefix <- args[3]

samplesheet <- read_csv(samplesheet) %>%
  # Get distinct samples, removing extra sequencing reps
  distinct(Sample, .keep_all = TRUE)

if (prefix != "ALL"){
  samplesheet <- samplesheet %>% filter(Condition == prefix)
}

CPMadjTMM_invNorm <- read_delim(normQuant, delim = "\t")

CPMadjTMM_invNorm_info <- CPMadjTMM_invNorm %>%
  dplyr::select(-`#chr`, -start, -end, -strand, -gene_name) %>%
  t() %>% as.data.frame() %>% row_to_names(row_number = 1)

if (prefix == "ALL"){
  CPMadjTMM_invNorm_info <- CPMadjTMM_invNorm_info %>% rownames_to_column("Sample") %>%
    left_join(samplesheet)
} else {
  CPMadjTMM_invNorm_info <- CPMadjTMM_invNorm_info %>% rownames_to_column("Donor") %>%
    left_join(samplesheet)
}

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

if (prefix != "ALL"){
  # Write all sample/PC's/sdev to file for correlation
  pc_sdev <- rbind(as.data.frame(CPMadjTMM_invNorm_pca$x), CPMadjTMM_invNorm_pca$sdev)
  rownames(pc_sdev)[length(rownames(pc_sdev))] <- "sdev"
  
  pc_sdev <- rownames_to_column(pc_sdev, var = "Donor")
  
  write_csv(pc_sdev, file = paste0("output/pca/", prefix, "_PC.csv"))
}

# Write sample/PC's to file for covariates
write_csv(bind_cols(CPMadjTMM_invNorm_info[,1], CPMadjTMM_invNorm_pca$x), file = paste0("output/covar/", prefix, "_PC.csv"))
