#!/usr/bin/R
library(readr)
library(dplyr)
library(janitor)
library(peer)

args <- commandArgs(trailingOnly = TRUE)
# Function to run PEER for different numbers of confounders

runPEER <- function(Nk, data, prefix){
  # Initialize model
  model <- PEER()
  PEER_setPhenoMean(model, data)
  
  # Optional: set observed covariates
  #PEER_setCovariates(model, as.matrix(covariates))
  
  PEER_setNmax_iterations(model, 1000)
  
  # Set number of confounders to infer
  PEER_setNk(model, Nk)
  
  # Run model
  # bound = 0.001, variance = 0.00001, Alpha a = 0.001, Alpha b = 0.1, Eps a = 0.1, Eps b = 10
  PEER_update(model)
  
  # Get output
  factors <- as.data.frame(PEER_getX(model))
  colnames(factors) <- paste0("PEER", seq(1, Nk))
  factors$Donor <- rownames(data)
  
  # Write to output
  write_csv(factors, file = paste0("output/covar/", prefix,"_PEERfactors_k", Nk, ".txt"))
  
  
  # Check variance of residuals to see if its decreasing and we want to add more iterations
  residualVariance <- as.data.frame(PEER_getResidualVars(model))
  write_csv(residualVariance, file = paste0("output/covar/", prefix, "_PEERresidualvariance_k", Nk, ".txt"))
  
  
}

# Read in data
expressionData <- read_delim(args[1]) %>%
  dplyr::select(-`#chr`, -start, -end, -strand, -gene_name) %>%
  t() %>% as.data.frame() %>% row_to_names(row_number = 1) %>%
  mutate(across(where(is.character), as.numeric)) %>%
  as.matrix()
  
lapply(seq(10, 100, 10), runPEER, data = expressionData, prefix = args[2])

