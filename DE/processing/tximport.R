#!/usr/bin/R
library(tximeta)
library(readr)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)

# IMPORT QUANTIFICATION WITH TXIMETA -------------------------------------------

# Read in sample sheet for coldata
coldata <- read_csv(args[1]) %>%
  # Get distinct samples
  distinct(Sample, .keep_all = TRUE) %>%
  ## Rename columns
  dplyr::rename(names = Sample)


# Add sample quants
coldata$files <- file.path("output", "quant", coldata$names, "quant.sf")

# Import salmon transcript quantification---------------------------------------
se <- tximeta(coldata)

# Convert to gene-level scaled transcripts -------------------------------------
gse <- summarizeToGene(se)

save(gse, file = paste0("output/", Sys.Date(), "_gse.rda"))