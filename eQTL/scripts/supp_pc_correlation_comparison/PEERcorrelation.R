library(tidyverse)
library(lubridate)
source("scripts/utils.R")
args <- commandArgs(trailingOnly = TRUE)
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

# Read in PEER factors and join with covariate information
if (condition != "ALL"){
  covariates <- covariates %>% filter(Condition == condition)
  peer <- read_csv(args[1]) %>%
    left_join(covariates, by = "Donor")
} else {
  peer <- read_csv(args[1]) %>%
    rename(Sample = Donor) %>%
    left_join(covariates, by = "Sample")
}

# Convert covariates to factor then numeric, and filter out any that only have
# 1 unique value
peer <- peer %>% 
  mutate(across(where(is.character), as.factor)) %>%
  mutate(across(where(is.factor), as.numeric)) %>%
  mutate(across(where(is.Date), as.numeric)) %>%
  mutate(across(where(is.POSIXct), as.numeric)) %>%
  select_if(~n_distinct(.) > 1) 

# Pearson's correlations and p-values
cors <- correlationTests_cor(peer, peer) %>% 
  as.data.frame() %>%
  rownames_to_column(var = "covariate")
  
cor_p <- correlationTests_p(peer, peer) %>% 
  as.data.frame() %>%
  rownames_to_column(var = "covariate")

# Write to file
write_csv(cors, file = paste0("output/peer/", condition, "_PEERcorrelations.csv"))
write_csv(cor_p, file = paste0("output/peer/", condition, "_PEERcorrelations_pval.csv"))

