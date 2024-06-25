library(tidyverse)
library(data.table)
library(plotgardener)
library(scales)
library(GenomicRanges)
library(GenomicFeatures)
library(mariner)
library(InteractionSet)
library(nullranges)
library(org.Hs.eg.db)
library(ggtext)
library(stats)
source("../plotting_utils.R")
source("../utils.R")

make_lead_eGene_gi <- function(eSNP_eGene, promoters_mapped_genes){
  eGene_id <- eSNP_eGene[["gene_id"]]
  chrom <- eSNP_eGene[["gene_chr"]]
  lead_var_pos <- as.numeric(eSNP_eGene[["variant_start"]])
  variantID <- eSNP_eGene[["variantID"]]
  print(variantID)
  
  # Make GRanges for lead variant
  lead_var_gr <- GRanges(seqnames = chrom,
                         ranges = IRanges(start = lead_var_pos,
                                          end = lead_var_pos))
  
  
  # Grab eGene_id promoter ranges
  eGene_promoters <- promoters_mapped_genes |> 
    filter(GENEID == eGene_id)
  
  if (length(eGene_promoters) > 0){
    # Create GInteractions between the lead and these promoters
    lead_eGene_gr <- GInteractions(anchor1 = rep(lead_var_gr, length(eGene_promoters)), 
                                   anchor2 = eGene_promoters,
                                   variantID = variantID)
    # Rename anchor2.GENEID column and anchor2.TXID column
    lead_eGene_gr$gene_id <- lead_eGene_gr$anchor2.GENEID
    lead_eGene_gr$anchor2.GENEID <- NULL
    lead_eGene_gr$tx_id <- lead_eGene_gr$anchor2.TXID
    lead_eGene_gr$anchor2.TXID <- NULL
  } else {
    # If no gene info, create empty GInteractions
    lead_eGene_gr <- GInteractions(variantID = variantID, 
                                   gene_id = eGene_id, tx_id = eGene_id)
  }
  
  
  
  return(lead_eGene_gr)
  
}

#### TXDB PARSING
# Get ENSEMBL Txdb
ensembl_txdb <- makeTxDbFromEnsembl(organism = "Homo sapiens", release = 111)
# Only keep standard chromosomes for compatibility with Hi-C files
ensembl_txdb <- keepStandardChromosomes(ensembl_txdb)
# Update seqlevels to UCSC style for compatibility with eqtl data
seqlevelsStyle(ensembl_txdb) <- "UCSC"
# Get ENSEMBL TxDb gene promoters
ensembl_promoter_genes <- promoters(ensembl_txdb, columns = c("GENEID", "TXID"))
# convert GENEID CharacterList just to character
ensembl_promoter_genes$GENEID <- as.character(ensembl_promoter_genes$GENEID)


# Shared eGenes
pbs_egenes <- read_csv("data/eqtl/CTL_PEER_k20_genoPC_perm1Mb_sig_rsID.csv") |> 
  pull(gene_id)
fnf_egenes <- read_csv("data/eqtl/FNF_PEER_k22_genoPC_perm1Mb_sig_rsID.csv") |> 
  pull(gene_id)

shared_egenes <- intersect(pbs_egenes, fnf_egenes)

shared_signals_eSNPs <- list()

for (chrom in 1:22){
  
  # Shared eqtls from either pbs or fnf
  shared_chrom_signal_eqtls <- bind_rows(fread(paste0("data/eqtl/qtl_cond/CTL_PEER_k20_genoPC_allSignals_nom1Mb_MAFs_chr",
                                                      chrom, ".csv"), data.table = FALSE) |> 
                                           filter(gene_id %in% shared_egenes) |> 
                                           filter(nom_sig == 1) |> 
                                           dplyr::select(gene_id, variantID, gene_chr, variant_start),
                                         fread(paste0("data/eqtl/qtl_cond/FNF_PEER_k22_genoPC_allSignals_nom1Mb_MAFs_chr",
                                                      chrom, ".csv"), data.table = FALSE) |> 
                                           filter(gene_id %in% shared_egenes) |>
                                           filter(nom_sig == 1) |> 
                                           dplyr::select(gene_id, variantID, gene_chr, variant_start)) |> 
    distinct()
  
  shared_signals_eSNPs[[paste0("chr", chrom)]] <- shared_chrom_signal_eqtls
  
}

shared_signals_eSNPs <- bind_rows(shared_signals_eSNPs)

shared_signals_eSNPs_eGenes <- apply(shared_signals_eSNPs, 1, make_lead_eGene_gi, ensembl_promoter_genes)
shared_signals_eSNPs_eGenes <- do.call(c, shared_signals_eSNPs_eGenes)

shared_signals_eSNPs_eGenes_10kb <- snapToBins(shared_signals_eSNPs_eGenes, 
                                               binSize = 10000)
save(shared_signals_eSNPs_eGenes_10kb, file = "data/hic/shared_signals_eSNPs_eGenes_10kb.rda")