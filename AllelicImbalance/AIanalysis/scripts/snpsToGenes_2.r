library(GenomicRanges)
library(data.table)
library(dplyr)
library(AnnotationDbi)
library(tools)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
source('scripts/utils.R')


# Parsing variantIDs ------------------------------------------------------
# Load DESeq results
## CTL
load('data/2023-01-10_AIresCTL.rda')
# Get GRanges of variant positions
AIresCTL_ranges <- varID_to_GRanges(rownames(notNA_resCTL))

## FNF
load('data/2023-01-10_AIresFNF.rda')
AIresFNF_ranges <- varID_to_GRanges(rownames(notNA_resFNF))

# Converting positions to genes -------------------------------------------
# Load TxDb and orgDb needed for gene conversions
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
orgdb <- org.Hs.eg.db

# CTL variants
AIresCTL_genes <- GRanges_to_Genes(ranges = AIresCTL_ranges, txdb = txdb, orgdb = orgdb)
write.csv(AIresCTL_genes, file = paste0("data/", Sys.Date(), "_AIresCTL_genes.csv"),
          quote = FALSE, row.names = FALSE)

# FNF variants
AIresFNF_genes <- GRanges_to_Genes(ranges = AIresFNF_ranges, txdb = txdb, orgdb = orgdb)
write.csv(AIresFNF_genes, file = paste0("data/", Sys.Date(), "_AIresFNF_genes.csv"),
          quote = FALSE, row.names = FALSE)
