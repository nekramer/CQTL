#!/usr/bin/R
library(tximeta)
library(readr)
library(dplyr)
library(tibble)
library(stringr)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(DESeq2)

# DIFFERENTIAL EXPRESSION ANALYSIS ---------------------------------------------

pval <- 0.01
log2FC <- 2


#load('output/')


## Convert to factors (avoids a warning)
colData(gse)[] <- lapply(colData(gse), factor)

## Build DESeq2 object
dds <- DESeqDataSet(gse, design = ~Donor + Condition)

## Filter out lowly expressed genes
# 10 reads in at least 5% of samples

keep <- rowSums(counts(dds) >= 10) >= ceiling(nrow(coldata)*0.05)
dds <- dds[keep,]

## Fit model
dds <- DESeq(dds)
save(dds, file = "data/differential_expression_dds.rda")
