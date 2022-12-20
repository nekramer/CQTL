#!/usr/bin/R
library(tidyverse)
library(lubridate)
source("scripts/utils.R")

args <- commandArgs(trailingOnly = TRUE)
condition <- args[2]
corOutput <- args[3]
pOutput <- args[4]

# Separate out pcs from covariates and select top 10
pcs <- read_csv(args[1]) %>%
  dplyr::select(contains("PC"))

# Read in covariates and convert to numerics
covariates <- read_csv(args[1]) %>%
  dplyr::select(!contains("PC")) %>%
  mutate(across(where(is.character), as.factor)) %>%
  mutate(across(where(is.factor), as.numeric)) %>%
  mutate(across(where(is.Date), as.numeric)) %>%
  mutate(across(where(is.POSIXct), as.numeric)) %>%
  # Make sure enough finite observations
  dplyr::select(where(~sum(is.na(.x)) < 3))
  
# Pearson's correlations and p-values
cors <- correlationTests_cor(pcs, covariates) %>% 
  na.omit() %>% 
  as.data.frame() %>%
  rownames_to_column(var = "covariate")
cor_p <- correlationTests_p(pcs, covariates) %>% 
  na.omit() %>%
  as.data.frame() %>%
  rownames_to_column(var = "covariate")

# Write to file
write_csv(cors, file = corOutput)
write_csv(cor_p, file = pOutput)