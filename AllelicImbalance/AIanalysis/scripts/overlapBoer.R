library(DESeq2)
library(tidyverse)
library(data.table)
library(prodlim)
source("scripts/utils.R")
# Read in significant AI variants and their LD buddies ------------------------

CTL_AI_sig <- fread("data/AIsigCTL05_final.csv", data.table = FALSE,
                    drop = c("count", "condition", "donor", "allele")) %>%
  distinct(.keep_all = TRUE)
FNF_AI_sig <- fread("data/AIsigFNF05_final.csv", data.table = FALSE,
                    drop = c("count", "condition", "donor", "allele")) %>%
  distinct(.keep_all = TRUE)


CTL_snps_only <- CTL_AI_sig[which(!CTL_AI_sig$rsid %in% FNF_AI_sig$rsid),]
FNF_snps_only <- FNF_AI_sig[which(!FNF_AI_sig$rsid %in% CTL_AI_sig$rsid),]



# # Compare to Boer nominal stats -------------------------------------------
# 
OAsubtypes <- list.files("/work/users/n/e/nekramer/External/gwas/Boer_2021_hg19")
OAsubtypes <- OAsubtypes[OAsubtypes != 'plots']
OAsubtypes <- OAsubtypes[OAsubtypes != 'Boer_effector_genes.csv']
OAsubtypes <- OAsubtypes[OAsubtypes != 'Boer_et_al_2021_S10.csv']


for (subtype in OAsubtypes){

  for (chr in 1:22){

    summary_stats <- fread(paste0("/work/users/n/e/nekramer/External/gwas/Boer_2021_hg19/",
                                  subtype, "/summary_statistics/", subtype,
                                  "_chr", chr, "_summary05Final.csv"),
                           data.table = FALSE,
                           drop = c("snp", "chrom", "pos")) %>%
      dplyr::rename("GWAS_rsid" = "rsid")

    # Filter for either overlapping AI variant or LD buddy of AI variant
    CTL_Boer <- CTL_snps_only %>%
      filter(rsid %in% summary_stats$GWAS_rsid | SNP_B_rsid %in% summary_stats$GWAS_rsid)
    FNF_Boer <- FNF_snps_only %>%
      filter(rsid %in% summary_stats$GWAS_rsid | SNP_B_rsid %in% summary_stats$GWAS_rsid)


    # Join back GWAS p
    CTL_Boer_gwasstats <- left_join(CTL_Boer, summary_stats,
                                    by = c("rsid" = "GWAS_rsid"), keep = TRUE) %>%
      left_join(summary_stats, by = c("SNP_B_rsid" = "GWAS_rsid"), keep = TRUE)
    FNF_Boer_gwasstats <- left_join(FNF_Boer, summary_stats,
                                    by = c("rsid" = "GWAS_rsid"), keep = TRUE) %>%
      left_join(summary_stats, by = c("SNP_B_rsid" = "GWAS_rsid"), keep = TRUE)
    # Append to files
    write_csv(CTL_Boer_gwasstats,
              file = paste0("data/", subtype, "_CTL_AI_LD_nominalOverlaps.csv"),
              append = TRUE)
    write_csv(FNF_Boer_gwasstats,
              file = paste0("data/", subtype, "_FNF_AI_LD_nominalOverlaps.csv"),
              append = TRUE)

  }
}

ctlOverlaps <- list()
fnfOverlaps <- list()

for (subtype in OAsubtypes){

  ctl <- read_csv(paste0("data/", subtype, "_CTL_AI_LD_nominalOverlaps.csv"),
                  col_names = c("rsid", "baseMean", "log2FoldChange", "lfcSE",
                                "stat", "pvalue", "padj", "chr", "pos",
                                "ref", "alt", "gene_symbol", "gene_start",
                                "gene_end", "gene_strand", "ENTREZID",
                                "variantID", "SNP_B", "R2", "SNP_B_rsid",
                                "EffectAllele.x", "AlternateAllele.x",
                                "EffectAlleleFrequency.x", "EffectSize.Beta.x",
                                "p.x", "SampleSize.x", "EffectSize.OR.x",
                                "GWAS_rsid.x", "EffectAllele.y", "AlternateAllele.y",
                                "EffectAlleleFrequency.y", "EffectSize.Beta.y",
                                "p.y", "SampleSize.y", "EffectSize.OR.y",
                                "GWAS_rsid.y"),
                  col_types = c("cddddddcdcccddccccdcccdddddcccdddddc")) %>%
    distinct() %>%
    mutate(OAsubtype = subtype)

  
  
  fnf <- read_csv(paste0("data/", subtype, "_FNF_AI_LD_nominalOverlaps.csv"),
                  col_names = c("rsid", "baseMean", "log2FoldChange", "lfcSE",
                                "stat", "pvalue", "padj", "chr", "pos",
                                "ref", "alt", "gene_symbol", "gene_start",
                                "gene_end", "gene_strand", "ENTREZID",
                                "variantID", "SNP_B", "R2", "SNP_B_rsid",
                                "EffectAllele.x", "AlternateAllele.x",
                                "EffectAlleleFrequency.x", "EffectSize.Beta.x",
                                "p.x", "SampleSize.x", "EffectSize.OR.x",
                                "GWAS_rsid.x", "EffectAllele.y", "AlternateAllele.y",
                                "EffectAlleleFrequency.y", "EffectSize.Beta.y",
                                "p.y", "SampleSize.y", "EffectSize.OR.y",
                                "GWAS_rsid.y"),
                  col_types = c("cddddddcdcccddccccdcccdddddcccdddddc")) %>%
    distinct() %>%
    mutate(OAsubtype = subtype)
  
  ctlOverlaps[[subtype]] <- ctl
  fnfOverlaps[[subtype]] <- fnf

}

ctlOverlaps <- bind_rows(ctlOverlaps) %>%
  arrange(padj)
write_csv(ctlOverlaps, file = "data/Boer_CTL_AI_LD_nominalOverlaps.csv")

fnfOverlaps <- bind_rows(fnfOverlaps) %>%
  arrange(padj)
write_csv(fnfOverlaps, file = "data/Boer_FNF_AI_LD_nominalOverlaps.csv")
