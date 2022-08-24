library(GenomicRanges)
library(data.table)
library(dplyr)
library(AnnotationDbi)
library(tools)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
source('scripts/utils.R')


# Parsing variantIDs ------------------------------------------------------

## CTL
# Load significant CTL variants
load('data/2022-08-23_AIsigCTL.rda')
# Get GRanges of variant positions
AIsigCTL_ranges <- varID_to_GRanges(rownames(sigCTL))

## FNF
load('data/2022-08-23_AIsigFNF.rda')
AIsigFNF_ranges <- varID_to_GRanges(rownames(sigFNF))

# Converting positions to genes -------------------------------------------
# Load TxDb and orgDb needed for gene conversions
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
orgdb <- org.Hs.eg.db

# CTL variants
AIsigCTL_genes <- GRanges_to_Genes(ranges = AIsigCTL_ranges, txdb = txdb, orgdb = orgdb)
write.csv(AIsigCTL_genes, file = paste0("data/", Sys.Date(), "_AIsigCTL_genes.csv"),
          quote = FALSE, row.names = FALSE)

# FNF variants
AIsigFNF_genes <- GRanges_to_Genes(ranges = AIsigFNF_ranges, txdb = txdb, orgdb = orgdb)
write.csv(AIsigFNF_genes, file = paste0("data/", Sys.Date(), "_AIsigFNF_genes.csv"),
          quote = FALSE, row.names = FALSE)
