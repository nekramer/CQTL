library(stringr)
library(purrr)
library(SNPlocs.Hsapiens.dbSNP155.GRCh38)
library(BSgenome)

source('scripts/utils.R')

# Parsing variantIDs ------------------------------------------------------

## CTL

# Load data
load('data/2022-08-17_AIresCTL.rda')

# Get GRanges of variant positions
AIresCTL_ranges <- varID_to_GRanges(rownames(notNA_resCTL))

## FNF



# Map positions to rsids --------------------------------------------------

dbSNP <- SNPlocs.Hsapiens.dbSNP155.GRCh38
