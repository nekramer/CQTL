library(tidyverse)
library(data.table)
library(plotgardener)
library(scales)
library(GenomicRanges)
library(mariner)
library(InteractionSet)
library(nullranges)
library(org.Hs.eg.db)
library(ggtext)
library(stats)
library(GenomicFeatures)
source("../plotting_utils.R")
source("../utils.R")
library(optparse)


args <- commandArgs(trailingOnly = TRUE)
chrom <- args[1]


determine_ld_loop_support <- function(eSNP_eGene, 
                                      promoters_mapped_genes, 
                                      loops,
                                      loop_binsize){
  eGene_id <- eSNP_eGene[["gene_id"]]
  variant_pos <- as.numeric(eSNP_eGene[["ld_pos"]])
  chrom <- eSNP_eGene[["variant_chr"]]
  
  # Make variant a GRanges object for overlapping
  eSNP_gr <- GRanges(seqnames = chrom,
                     ranges = IRanges(start = variant_pos,
                                      end = variant_pos))
  
  # Grab eGene_id promoter ranges
  eGene_promoters <- promoters_mapped_genes |> 
    filter(GENEID == eGene_id)
  
  # Get loop interactions that link snp and any of the eGene's promoters
  loop_interactions <- linkOverlaps(snapToBins(loops, binSize = loop_binsize), 
                                    subject1 = eSNP_gr, 
                                    subject2 = eGene_promoters)
  
  if (nrow(loop_interactions) > 0){
    return(names(loops[unique(loop_interactions$query)]))
  } else {
    return("no_loop")
  }
  
}

# PBS eQTLs on chrom and LD > 0.8
pbs_eqtls_LD_chrom <- fread("data/eqtl/CTL_PEER_k20_genoPC_cond1Mb_topSignals_rsID_LD.csv",
                      data.table = FALSE) |> 
  filter(variant_chr == chrom) |> 
  filter(R2 > 0.6) |> 
  separate_wider_delim(cols = "ld_variantID", delim = ":", 
                       names = c(NA, "ld_pos", NA, NA), 
                       cols_remove = FALSE)

# FNF eQTLs on chrom and LD > 0.8
fnf_eqtls_LD_chrom <- fread("data/eqtl/FNF_PEER_k22_genoPC_cond1Mb_topSignals_rsID_LD.csv",
                      data.table = FALSE) |> 
  filter(variant_chr == chrom) |> 
  filter(R2 > 0.6) |> 
  separate_wider_delim(cols = "ld_variantID", delim = ":", 
                       names = c(NA, "ld_pos", NA, NA), 
                       cols_remove = FALSE)

# Bind PBS and FN-f
PBS_FNF_eqtls_LD_chrom <- bind_rows(pbs_eqtls_LD_chrom, fnf_eqtls_LD_chrom) |> 
  distinct(gene_id, ld_variantID, .keep_all = TRUE)

# Load static and sig gained loops and concat
load("data/hic/CQTL_10kb_static_loops.rda")
load("data/hic/CQTL_10kb_sig_gained_loops.rda")
CQTL_10kb_allLoops <- c(CQTL_10kb_static_loops, CQTL_10kb_sig_gained_loops)
# Add names to loops for easier access
names(CQTL_10kb_allLoops) <- c(paste0("static_loop_", 
                                      seq(1:length(CQTL_10kb_static_loops))),
                               paste0("gained_loop_", 
                                      seq(1:length(CQTL_10kb_sig_gained_loops))))

# Filter for chrom loops
CQTL_10kb_allLoops_chrom <- CQTL_10kb_allLoops[seqnames1(CQTL_10kb_allLoops) == chrom]


# Get ENSEMBL TxDb gene promoters
ensembl_txdb <- makeTxDbFromEnsembl(organism = "Homo sapiens", release = 111)
# Update seqlevels to UCSC style for compatibility with eqtl data
seqlevelsStyle(ensembl_txdb) <- "UCSC"
ensembl_promoter_genes <- promoters(ensembl_txdb, columns = c("GENEID", "TXID"))
# convert GENEID CharacterList just to character
ensembl_promoter_genes$GENEID <- as.character(ensembl_promoter_genes$GENEID)

# Add additional column for loop support
PBS_FNF_eqtls_LD_chrom$looping <- apply(PBS_FNF_eqtls_LD_chrom, 
                                  1, 
                                  determine_ld_loop_support, 
                                  ensembl_promoter_genes,
                                  CQTL_10kb_allLoops_chrom,
                                  loop_binsize = 10000)

write_csv(PBS_FNF_eqtls_LD_chrom, file = paste0("data/hic/PBS_FNF_eQTLs_LD_", chrom,".csv"))

