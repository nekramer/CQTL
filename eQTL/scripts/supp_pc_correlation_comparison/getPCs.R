#!/usr/bin/R
library(tidyverse)
library(janitor)

args <- commandArgs(trailingOnly = TRUE)

normQuant <- args[1]
condition <- args[2]
samplesheet <- args[3]
dnaSamplesheet <- args[4]
donorSamplesheet <- args[5]

# Read in samplesheets and join for covariate information
covariates <- read_csv(samplesheet) %>%
  # Get distinct samples, removing extra sequencing reps
  distinct(Sample, .keep_all = TRUE) %>%
  # Join with dnaSamplesheet by Donor
  left_join(read_csv(dnaSamplesheet) %>% dplyr::select(-Sample)) %>%
  # Join with donorSamplesheet by Donor
  left_join(read_csv(donorSamplesheet)) %>%
  # Remove any columns where all are NA 
  discard(~all(is.na(.x)))

if (condition != "ALL"){
  covariates <- covariates %>% filter(Condition == condition)
}



# Read in normalized counts
CPMadjTMM_invNorm <- read_delim(normQuant, delim = "\t") %>%
  column_to_rownames(var = "gene_id") %>%
  # Filter out gene_info
  dplyr::select(-`#chr`, -start, -end, -strand, -gene_name)

# Calculate variances to remove columns where variance is 0 (can't do PCA)
variances <- CPMadjTMM_invNorm %>% 
  t() %>% apply(2, var, na.rm = TRUE)

# Get PC's
if (length(which(variances == 0) > 0)){
  
  CPMadjTMM_invNorm_pca <- CPMadjTMM_invNorm %>%
    t() %>% as.data.frame() %>% dplyr::select(-which(variances == 0)) %>%
    prcomp(center = TRUE, scale = TRUE)
  
} else {
  
  CPMadjTMM_invNorm_pca <- CPMadjTMM_invNorm %>%
    t() %>%
    prcomp(center = TRUE, scale = TRUE)
}


# Variance Explained ------------------------------------------------------

var_explained <- data.frame("PC" = paste0("PC", 1:length(CPMadjTMM_invNorm_pca$sdev)),
                            "sdev" = CPMadjTMM_invNorm_pca$sdev) %>%
  mutate(varianceExplained = sdev^2/sum(sdev^2)) %>%
  # Percent variance explained
  mutate(percentVarianceExplained = 100*varianceExplained) 

# Write to file
write_csv(var_explained, 
          file = paste0("output/pca/", condition, "_PC_varExplained.csv"))


# Final PCs ---------------------------------------------------------------

# Join PCs with covariate information
if (condition == "ALL"){
  final_pca <- as.data.frame(CPMadjTMM_invNorm_pca$x) %>% 
    rownames_to_column(var = "Sample") %>%
    left_join(covariates, by = "Sample")
  
} else {
  final_pca <- as.data.frame(CPMadjTMM_invNorm_pca$x) %>% 
    rownames_to_column(var = "Donor") %>%
    left_join(covariates, by = "Donor")
}

# Write to file
write_csv(final_pca, 
          file = paste0("output/pca/", condition, "_PCs.csv"))

