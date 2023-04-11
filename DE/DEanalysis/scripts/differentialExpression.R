#!/usr/bin/R
library(tximeta)
library(readr)
library(dplyr)
library(tibble)
library(stringr)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(DESeq2)

# DIFFERENTIAL EXPRESSION ANALYSIS ---------------------------------------------

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

## Convert results to GRanges
de_genes <- results(dds, name = "Condition_FNF_vs_CTL",
                    format = "GRanges") %>%
  plyranges::names_to_column("gene_id")

## Add gene symbols to de_genes and
## remove non-standard chromosomes
de_genes <- 
  inner_join(x = as.data.table(de_genes),
             y = as.data.table(rowData(gse)) %>%
               dplyr::select(c("gene_id", "symbol", "tx_ids")),
             by = "gene_id") %>%
  makeGRangesFromDataFrame(keep.extra.columns = TRUE) %>%
  keepStandardChromosomes(pruning.mode = "coarse")

## Add seqinfo to object
txdb <- 
  TxDb.Hsapiens.UCSC.hg38.knownGene %>%
  keepStandardChromosomes()

seqlevels(de_genes) <- seqlevels(txdb)
seqinfo(de_genes) <- seqinfo(txdb)

# Save object
saveRDS(de_genes, file = file.path("data", "fnfControl_diffgenes.rds"))
