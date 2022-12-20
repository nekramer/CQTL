#!/usr/bin/R
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
correction <- args[3]
threshold <- as.numeric(args[4])



CTL_sig_eGenes <- read_csv(args[1]) %>%
  filter(.[[correction]] <= threshold) %>%
  dplyr::select(gene_id, gene_name, gene_chr, gene_start, gene_end, gene_strand,
                variantID, variant_chr, variant_start, variant_end, beta,
                contains(correction))

FNF_sig_eGenes <- read_csv(args[2]) %>%
  filter(.[[correction]] <= threshold) %>%
  dplyr::select(gene_id, gene_name, gene_chr, gene_start, gene_end, gene_strand,
                variantID, variant_chr, variant_start, variant_end, beta,
                contains(correction))


CTL_only <- CTL_sig_eGenes[!CTL_sig_eGenes$gene_id %in% FNF_sig_eGenes$gene_id,]
write_csv(CTL_only, file = "output/reQTL/CTLonly_sig_eGenes.csv")

FNF_only <- FNF_sig_eGenes[!FNF_sig_eGenes$gene_id %in% CTL_sig_eGenes$gene_id,]
write_csv(FNF_only, file = "output/reQTL/FNFonly_sig_eGenes.csv")
