#!/usr/bin/R
library(tidyverse)
library(tools)

args <- commandArgs(trailingOnly = TRUE)
threshold <- as.numeric(args[1])
outputFile <- args[2]
eGene_stats <- list()

for (file in args[3:length(args)]){
  # Read in file
  data <- read_csv(file)
  
  if (nrow(data) > 0){
    
    # Filter for significance based on qvalue and count number of eGenes
    qvalSig <- data %>% dplyr::filter(qval <= threshold) %>% nrow()
    
    # Filter for significance based on FDR and count number of eGenes
    fdrSig <- data %>% dplyr::filter(FDR <= threshold) %>% nrow()
    
    fileName <- file_path_sans_ext(basename(file))
    fileName <- gsub("_perm1Mb", "", fileName)
    
    # Add to list
    eGene_stats[[fileName]] <- data.frame("qval" = c(qvalSig),
                                          "FDR" = c(fdrSig))
  } else {
    eGene_stats[[fileName]] <- data.frame("qval" = c(0),
                                          "FDR" = c(0))
  }
  

}

# Combine list into dataframe and write to file
bind_rows(eGene_stats, .id = "Covariate_Options") %>% 
  write_csv(outputFile)
