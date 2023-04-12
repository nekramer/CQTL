library(ggplot2)
library(DESeq2)
source("scripts/utils.R")
source("scripts/plotAIdata.R")
# Read in and format CTL/FNF AI rsids/genes ------------------------------------

CTL <- reformatRSIDdata("data/2023-01-10_AIsigCTL05_rsids.csv")
FNF <- reformatRSIDdata("data/2023-01-10_AIsigFNF05_rsids.csv")

# DESeq results --------------------------------------------------------
load("data/2023-01-10_AIsigCTL_.05.rda")
load("data/2023-01-10_AIsigFNF_.05.rda")
load("data/dds.rda")

# Reformat and join with rsid/gene info
sigCTL_df <- reformatDESeqdata(sigCTL_.05, CTL)
sigFNF_df <- reformatDESeqdata(sigFNF_.05, FNF)


# CTL/FNF/both comparisons -----------------------------------------------------

CTL_snps_only <- sigCTL_df[which(!sigCTL_df$rsid %in% sigFNF_df$rsid),]
FNF_snps_only <- sigFNF_df[which(!sigFNF_df$rsid %in% sigCTL_df$rsid),]
both_snps_only <- sigCTL_df[which(sigCTL_df$rsid %in% sigFNF_df$rsid),]


# Comparison with differential genes --------------------------------------

sig_de_genes <- read_csv("data/sig_differential_genes.csv")

CTL_snps_deoverlap <- CTL_snps_only[which(CTL_snps_only$gene_symbol %in% sig_de_genes$symbol),]
FNF_snps_deoverlap <- FNF_snps_only[which(FNF_snps_only$gene_symbol %in% sig_de_genes$symbol),]





both_snps_deoverlap <- both_snps_only[which(both_snps_only$gene_symbol %in% sig_de_genes$symbol),]
