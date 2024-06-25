library(tidyverse)
library(vcfR)
source("../utils.R")
source("../plotting_utils.R")

# coloc'ed eGene boxplots - 
# TGFA, PBS-specific eQTL (PIK3R1), and FN-f-specific qtl (RNF144B)-----------------

pbs_egenes <- read_csv("data/eqtl/CTL_PEER_k20_genoPC_cond1Mb_topSignals_rsID.csv",
                       col_select = c("rsID", "gene_id", "variantID", "gene_symbol", "signal")) |>
  # all in 0 signal
  filter(gene_symbol %in% c("TGFA", "PIK3R1") & signal == 0) |> 
  mutate(grouping = ifelse(gene_symbol == "TGFA", "shared", "PBS-specific")) |> 
  mutate(source = "PBS")
fnf_egenes <- read_csv("data/eqtl/FNF_PEER_k22_genoPC_cond1Mb_topSignals_rsID.csv",
                       col_select = c("rsID", "gene_id", "variantID", "gene_symbol", "signal")) |> 
  filter(gene_symbol %in%  c("TGFA", "RNF144B") & signal == 0) |> 
  mutate(grouping = ifelse(gene_symbol == "TGFA", "shared", "FN-f-specific")) |> 
  mutate(source = "FNF")

boxplot_egenes <- bind_rows(pbs_egenes, fnf_egenes)

# Get minor allele for correct genotype ordering
boxplot_egenes$MA <- unlist(lapply(boxplot_egenes$variantID, get_minor_allele, boxplot_egenes))

# Get PBS and FN-f betas and nominal p-values
boxplot_egenes <- boxplot_egenes |> 
  bind_cols(lapply(boxplot_egenes$variantID, get_pvals_betas, boxplot_egenes) |> 
              bind_rows())

# Geno data from VCF file containing FNF lead eQTL variants
vcf <- vcfR2tidy(read.vcfR("data/vcf/CTLk20_FNFk22_genoPC_leadVars.vcf.gz", 
                           verbose = FALSE))
variants <- vcf$fix |> 
  dplyr::select(POS, ID, REF, ALT)
geno_data <- vcf$gt |> 
  left_join(variants, by = "POS") |> 
  filter(ID %in% boxplot_egenes$variantID) |> 
  dplyr::select(Indiv, gt_GT_alleles, ID) |> 
  dplyr::rename(Donor = Indiv,
                variantID = ID)

# Normalized expression data
CTL_normQuant <- read_delim("data/rna/CTL_CPMadjTMM_invNorm.bed.gz") |> 
  filter(gene_id %in% boxplot_egenes$gene_id) |> 
  dplyr::select(-`#chr`, -start, -end, -length, -strand) |> 
  pivot_longer(cols = starts_with("AM"), names_to = "Donor", values_to = "expression") |> 
  mutate(Condition = "PBS")

FNF_normQuant <- read_delim("data/rna/FNF_CPMadjTMM_invNorm.bed.gz") |> 
  filter(gene_id %in% boxplot_egenes$gene_id)  |> 
  dplyr::select(-`#chr`, -start, -end, -length, -strand) |> 
  pivot_longer(cols = starts_with("AM"), names_to = "Donor", values_to = "expression") |> 
  mutate(Condition = "FN-f")

ALL_normQuant <- bind_rows(CTL_normQuant, FNF_normQuant)

# Join all data elements together
egene_boxplot_data <- left_join(boxplot_egenes, geno_data, by = "variantID") |> 
  left_join(ALL_normQuant, by = c("gene_id", "Donor"), relationship = "many-to-many")

# Get orders of genotypes  
geno_orders <- lapply(unique(egene_boxplot_data$variantID), 
                      determine_geno_order, 
                      data = egene_boxplot_data) |> 
  bind_rows()

egene_boxplot_data <- egene_boxplot_data |> 
  left_join(geno_orders, by = c("variantID", "gt_GT_alleles"))
egene_boxplot_data$Condition <- factor(egene_boxplot_data$Condition, levels = c("PBS", "FN-f"))

# Create boxplots
coloc_shared_eqtl_boxplot <- create_eqtl_boxplot_horizontal(egene_boxplot_data,
                                                            group = "shared",
                                                            stat_loc = "top_right",
                                                            rsID_loc = "bottom")
coloc_pbs_eqtl_boxplot <- create_eqtl_boxplot_horizontal(egene_boxplot_data, 
                                                         group = "PBS-specific",
                                                         stat_loc = "top_right",
                                                         rsID_loc = "bottom")
coloc_fnf_eqtl_boxplot <- create_eqtl_boxplot_horizontal(egene_boxplot_data, 
                                                         group = "FN-f-specific",
                                                         stat_loc = "top_right",
                                                         rsID_loc = "bottom")

save(coloc_shared_eqtl_boxplot, 
     file = "plots/colocs_Fig5/coloc_shared_eqtl_boxplot.rda")
save(coloc_pbs_eqtl_boxplot, 
     file = "plots/colocs_Fig5/coloc_pbs_eqtl_boxplot.rda")
save(coloc_fnf_eqtl_boxplot, 
     file = "plots/colocs_Fig5/coloc_fnf_eqtl_boxplot.rda")

# Data for locus zooms ----------------------------------------------------
## TGFA
tgfa_egene <- read_csv("data/eqtl/CTL_PEER_k20_genoPC_cond1Mb_topSignals_rsID.csv") |>
  filter(gene_symbol == "TGFA" & signal == 0)
tgfa_egene_ld <- fread("data/eqtl/CTL_PEER_k20_genoPC_cond1Mb_topSignals_rsID_LD.csv",
                      data.table = FALSE) |> 
  filter(gene_symbol == "TGFA" & signal == 0) |> 
  dplyr::select(ld_variantID, R2)

# GWAS
tgfa_coloc <- read_csv("data/colocalization/processed/PBS_FNF_allSignals_EUR_GWAS_250Kb_colocs.csv") |> 
  filter(eGene_name == "TGFA")
tgfa_qtl_rsid <- tgfa_coloc |> pull(qtl_rsid)
tgfa_gwas_rsid <- tgfa_coloc |> pull(GWAS_lead)

# Lead information
gwas_tgfa_info <- read_csv(paste0("/work/users/n/e/nekramer/External/gwas/Boer_reprocessed/",
                                 tgfa_coloc$OA, "/leads/EUR_", 
                                 tgfa_coloc$OA, "_leads_ld_final.csv"),
                          col_types = "cdddccccddddddddddddddddddddddddddddddddddddddddddccdcccd") |> 
  filter(ldbuddy_rsID == tgfa_qtl_rsid) |> 
  mutate(id_v1 = paste0("chr", `CHR:hg38POS`, ":", EA, ":", NEA),
         id_v2 = paste0("chr", `CHR:hg38POS`, ":", NEA, ":", EA))
tgfa_gwas_lead_ids <- unique(c(gwas_tgfa_info[["id_v1"]], gwas_tgfa_info[["id_v2"]]))

# Signal with LD buddies
gwas_tgfa_data <- read_csv(paste0("/work/users/n/e/nekramer/External/gwas/Boer_reprocessed/",
                                   tgfa_coloc$OA, "/leads/EUR_",
                                   tgfa_coloc$OA, "_leads_ld_final.csv"),
                           col_types = "cdddccccddddddddddddddddddddddddddddddddddddddddddccdcccd") |> 
  filter(rsID == tgfa_gwas_rsid) |> 
  filter(chrom == gsub("chr", "", tgfa_egene[["gene_chr"]])) |> 
  mutate(ldbuddy_hg38pos = as.integer(str_extract(`ldbuddy_CHR:hg38POS`, "(?<=:)[0-9]+"))) |> 
  filter(ldbuddy_hg38pos >= (as.numeric(tgfa_egene[["variant_start"]]) - 500000) 
         & ldbuddy_hg38pos <= (as.numeric(tgfa_egene[["variant_start"]]) + 500000)) |> 
  dplyr::select(ldbuddy_rsID, `ldbuddy_CHR:hg38POS`, chrom, ldbuddy_hg38pos, p, ldbuddy_R2) |> 
  dplyr::rename(snp = ldbuddy_rsID, pos = ldbuddy_hg38pos)

# Get summary stats for other variants in region w/o R2 info
tgfa_summary_stats <- fread(paste0("/work/users/n/e/nekramer/External/gwas/Boer_reprocessed/",
                                    tgfa_coloc$OA, "/summary_stats/", 
                                    tgfa_coloc$OA, "_", 
                                    tgfa_egene[["gene_chr"]], ".csv"),
                       data.table = FALSE) |> 
  filter(hg38pos >= (as.numeric(tgfa_egene[["variant_start"]]) - 500000) & 
           hg38pos <= (as.numeric(tgfa_egene[["variant_start"]]) + 500000)) |>
  filter(!`CHR:hg38POS` %in% gwas_tgfa_data$`ldbuddy_CHR:hg38POS`) |> 
  dplyr::select(`CHR:hg38POS`, chrom, hg38pos, p) |> 
  dplyr::rename(snp = `CHR:hg38POS`, pos = hg38pos) |> 
  mutate(ldbuddy_R2 = NA)

tgfa_gwas_signal <- gwas_tgfa_data |> 
  dplyr::select(-`ldbuddy_CHR:hg38POS`) |> 
  bind_rows(tgfa_summary_stats) |> 
  mutate(chrom = paste0("chr", chrom)) |> 
  # Create group for LD group
  mutate(LDgrp = factor(cut(ldbuddy_R2, c(0, 0.2, 0.4, 0.6, 0.8, 1)),
                        levels = c("(0,0.2]", "(0.2,0.4]", "(0.4,0.6]",
                                   "(0.6,0.8]", "(0.8,1]", NA)))
tgfa_gwas_signal$LDgrp <- addNA(tgfa_gwas_signal$LDgrp)

# eQTL
tgfa_signals <- bind_rows(fread(paste0("data/eqtl/qtl_nom/CTL_PEER_k20_genoPC_allSignals_nom1Mb_MAFs_", 
                                        tgfa_egene[["gene_chr"]], ".csv"),
                                 data.table = FALSE) |> 
                             filter(gene_id == tgfa_egene[["gene_id"]] & signal == tgfa_egene[["signal"]]) |> 
                             mutate(Condition = "PBS"),
                           fread(paste0("data/eqtl/qtl_nom/FNF_PEER_k22_genoPC_allSignals_nom1Mb_MAFs_", 
                                        tgfa_egene[["gene_chr"]], ".csv"),
                                 data.table = FALSE) |> 
                             filter(gene_id == tgfa_egene[["gene_id"]] & signal == tgfa_egene[["signal"]]) |> 
                             mutate(Condition = "FN-f")) |>
  mutate(rsID = ifelse(variantID == tgfa_egene[["variantID"]], tgfa_qtl_rsid,
                       ifelse(variantID %in% tgfa_gwas_lead_ids, tgfa_gwas_rsid, NA))) |> 
  dplyr::rename(chrom = variant_chr,
                pos = variant_start,
                p = nom_pval,
                snp = rsID) |> 
  left_join(tgfa_egene_ld, by = join_by(variantID == ld_variantID)) |> 
  mutate(LDgrp = cut(R2, c(0, 0.2, 0.4, 0.6, 0.8, 1)))
tgfa_signals$LDgrp <- addNA(tgfa_signals$LDgrp)

tgfa_pbs_qval <- read_csv("data/eqtl/qtl_supp/CTL_PEER_k20_genoPC_perm1Mb_FDR.csv") |> 
  filter(gene_symbol == "TGFA") |> 
  pull(qval)

tgfa_fnf_qval <- read_csv("data/eqtl/qtl_supp/FNF_PEER_k22_genoPC_perm1Mb_FDR.csv") |> 
  filter(gene_symbol == "TGFA") |> 
  pull(qval)

## PIK3R1
pik3r1_egene <- read_csv("data/eqtl/CTL_PEER_k20_genoPC_cond1Mb_topSignals_rsID.csv") |>
  filter(gene_symbol == "PIK3R1" & signal == 0)
pik3r1_egene_ld <- fread("data/eqtl/CTL_PEER_k20_genoPC_cond1Mb_topSignals_rsID_LD.csv",
                        data.table = FALSE) |> 
  filter(gene_symbol == "PIK3R1" & signal == 0) |> 
  dplyr::select(ld_variantID, R2)

# GWAS
pik3r1_coloc <- read_csv("data/colocalization/processed/PBS_FNF_allSignals_EUR_GWAS_250Kb_colocs.csv") |> 
  filter(eGene_name == "PIK3R1")
pik3r1_qtl_rsid <- pik3r1_coloc |> pull(qtl_rsid)
pik3r1_gwas_rsid <- pik3r1_coloc |> pull(GWAS_lead)

# Lead information
gwas_pik3r1_info <- read_csv(paste0("/work/users/n/e/nekramer/External/gwas/Boer_reprocessed/",
                                 pik3r1_coloc$OA, "/leads/EUR_", 
                                 pik3r1_coloc$OA, "_leads_ld_final.csv"),
                          col_types = "cdddccccddddddddddddddddddddddddddddddddddddddddddccdcccd") |> 
  filter(ldbuddy_rsID == pik3r1_qtl_rsid) |> 
  mutate(id_v1 = paste0("chr", `CHR:hg38POS`, ":", EA, ":", NEA),
         id_v2 = paste0("chr", `CHR:hg38POS`, ":", NEA, ":", EA))
pik3r1_gwas_lead_ids <- unique(c(gwas_pik3r1_info[["id_v1"]], gwas_pik3r1_info[["id_v2"]]))

# Signal with LD buddies
gwas_pik3r1_data <- read_csv(paste0("/work/users/n/e/nekramer/External/gwas/Boer_reprocessed/",
                                   pik3r1_coloc$OA, "/leads/EUR_",
                                   pik3r1_coloc$OA, "_leads_ld_final.csv"),
                           col_types = "cdddccccddddddddddddddddddddddddddddddddddddddddddccdcccd") |> 
  filter(rsID == pik3r1_gwas_rsid) |> 
  filter(chrom == gsub("chr", "", pik3r1_egene[["gene_chr"]])) |> 
  mutate(ldbuddy_hg38pos = as.integer(str_extract(`ldbuddy_CHR:hg38POS`, "(?<=:)[0-9]+"))) |> 
  filter(ldbuddy_hg38pos >= (as.numeric(pik3r1_egene[["variant_start"]]) - 500000) 
         & ldbuddy_hg38pos <= (as.numeric(pik3r1_egene[["variant_start"]]) + 500000)) |> 
  dplyr::select(ldbuddy_rsID, `ldbuddy_CHR:hg38POS`, chrom, ldbuddy_hg38pos, p, ldbuddy_R2) |> 
  dplyr::rename(snp = ldbuddy_rsID, pos = ldbuddy_hg38pos)

# Get summary stats for other variants in region w/o R2 info
pik3r1_summary_stats <- fread(paste0("/work/users/n/e/nekramer/External/gwas/Boer_reprocessed/",
                                    pik3r1_coloc$OA, "/summary_stats/", 
                                    pik3r1_coloc$OA, "_", 
                                    pik3r1_egene[["gene_chr"]], ".csv"),
                       data.table = FALSE) |> 
  filter(hg38pos >= (as.numeric(pik3r1_egene[["variant_start"]]) - 500000) & 
           hg38pos <= (as.numeric(pik3r1_egene[["variant_start"]]) + 500000)) |>
  filter(!`CHR:hg38POS` %in% gwas_pik3r1_data$`ldbuddy_CHR:hg38POS`) |> 
  dplyr::select(`CHR:hg38POS`, chrom, hg38pos, p) |> 
  dplyr::rename(snp = `CHR:hg38POS`, pos = hg38pos) |> 
  mutate(ldbuddy_R2 = NA)

pik3r1_gwas_signal <- gwas_pik3r1_data |> 
  dplyr::select(-`ldbuddy_CHR:hg38POS`) |> 
  bind_rows(pik3r1_summary_stats) |> 
  mutate(chrom = paste0("chr", chrom)) |> 
  # Create group for LD group
  mutate(LDgrp = factor(cut(ldbuddy_R2, c(0, 0.2, 0.4, 0.6, 0.8, 1)),
                        levels = c("(0,0.2]", "(0.2,0.4]", "(0.4,0.6]",
                                   "(0.6,0.8]", "(0.8,1]", NA)))
pik3r1_gwas_signal$LDgrp <- addNA(pik3r1_gwas_signal$LDgrp)

# eQTL
pik3r1_signals <- bind_rows(fread(paste0("data/eqtl/qtl_nom/CTL_PEER_k20_genoPC_allSignals_nom1Mb_MAFs_", 
                                         pik3r1_egene[["gene_chr"]], ".csv"),
                                  data.table = FALSE) |> 
                              filter(gene_id == pik3r1_egene[["gene_id"]] & signal == pik3r1_egene[["signal"]]) |> 
                              mutate(Condition = "PBS"),
                            fread(paste0("data/eqtl/qtl_nom/FNF_PEER_k22_genoPC_allSignals_nom1Mb_MAFs_", 
                                         pik3r1_egene[["gene_chr"]], ".csv"),
                                  data.table = FALSE) |> 
                              filter(gene_id == pik3r1_egene[["gene_id"]] & signal == pik3r1_egene[["signal"]]) |> 
                              mutate(Condition = "FN-f")) |>
  mutate(rsID = ifelse(variantID == pik3r1_egene[["variantID"]], pik3r1_qtl_rsid, 
                       ifelse(variantID %in% pik3r1_gwas_lead_ids, pik3r1_gwas_rsid, NA))) |> 
  dplyr::rename(chrom = variant_chr,
                pos = variant_start,
                p = nom_pval,
                snp = rsID) |> 
  left_join(pik3r1_egene_ld, by = join_by(variantID == ld_variantID)) |> 
  mutate(LDgrp = cut(R2, c(0, 0.2, 0.4, 0.6, 0.8, 1)))
pik3r1_signals$LDgrp <- addNA(pik3r1_signals$LDgrp)

pik3r1_pbs_qval <- read_csv("data/eqtl/qtl_supp/CTL_PEER_k20_genoPC_perm1Mb_FDR.csv") |> 
  filter(gene_symbol == "PIK3R1") |> 
  pull(qval)

pik3r1_fnf_qval <- read_csv("data/eqtl/qtl_supp/FNF_PEER_k22_genoPC_perm1Mb_FDR.csv") |> 
  filter(gene_symbol == "PIK3R1") |> 
  pull(qval)

## RNF144B
rnf144b_egene <- read_csv("data/eqtl/FNF_PEER_k22_genoPC_cond1Mb_topSignals_rsID.csv") |>
  filter(gene_symbol == "RNF144B" & signal == 0)
rnf144b_egene_ld <- fread("data/eqtl/FNF_PEER_k22_genoPC_cond1Mb_topSignals_rsID_LD.csv",
                         data.table = FALSE) |> 
  filter(gene_symbol == "RNF144B" & signal == 0) |> 
  dplyr::select(ld_variantID, R2)

# GWAS
rnf144b_coloc <- read_csv("data/colocalization/processed/PBS_FNF_allSignals_EUR_GWAS_250Kb_colocs.csv") |> 
  filter(eGene_name == "RNF144B")
rnf144b_qtl_rsid <- rnf144b_coloc |> pull(qtl_rsid)
rnf144b_gwas_rsid <- rnf144b_coloc |> pull(GWAS_lead)

# Lead information
gwas_rnf144b_info <- read_csv(paste0("/work/users/n/e/nekramer/External/gwas/Boer_reprocessed/",
                                    rnf144b_coloc$OA, "/leads/EUR_", 
                                    rnf144b_coloc$OA, "_leads_ld_final.csv"),
                             col_types = "cdddccccddddddddddddddddddddddddddddddddddddddddddccdcccd") |> 
  filter(ldbuddy_rsID == rnf144b_qtl_rsid) |> 
  mutate(id_v1 = paste0("chr", `CHR:hg38POS`, ":", EA, ":", NEA),
         id_v2 = paste0("chr", `CHR:hg38POS`, ":", NEA, ":", EA))
rnf144b_gwas_lead_ids <- unique(c(gwas_rnf144b_info[["id_v1"]], gwas_rnf144b_info[["id_v2"]]))

# Signal with LD buddies
gwas_rnf144b_data <- read_csv(paste0("/work/users/n/e/nekramer/External/gwas/Boer_reprocessed/",
                                    rnf144b_coloc$OA, "/leads/EUR_",
                                    rnf144b_coloc$OA, "_leads_ld_final.csv"),
                             col_types = "cdddccccddddddddddddddddddddddddddddddddddddddddddccdcccd") |> 
  filter(rsID == rnf144b_gwas_rsid) |> 
  filter(chrom == gsub("chr", "", rnf144b_egene[["gene_chr"]])) |> 
  mutate(ldbuddy_hg38pos = as.integer(str_extract(`ldbuddy_CHR:hg38POS`, "(?<=:)[0-9]+"))) |> 
  filter(ldbuddy_hg38pos >= (as.numeric(rnf144b_egene[["variant_start"]]) - 500000) 
         & ldbuddy_hg38pos <= (as.numeric(rnf144b_egene[["variant_start"]]) + 500000)) |> 
  dplyr::select(ldbuddy_rsID, `ldbuddy_CHR:hg38POS`, chrom, ldbuddy_hg38pos, p, ldbuddy_R2) |> 
  dplyr::rename(snp = ldbuddy_rsID, pos = ldbuddy_hg38pos)

# Get summary stats for other variants in region w/o R2 info
rnf144b_summary_stats <- fread(paste0("/work/users/n/e/nekramer/External/gwas/Boer_reprocessed/",
                                     rnf144b_coloc$OA, "/summary_stats/", 
                                     rnf144b_coloc$OA, "_", 
                                     rnf144b_egene[["gene_chr"]], ".csv"),
                              data.table = FALSE) |> 
  filter(hg38pos >= (as.numeric(rnf144b_egene[["variant_start"]]) - 500000) & 
           hg38pos <= (as.numeric(rnf144b_egene[["variant_start"]]) + 500000)) |>
  filter(!`CHR:hg38POS` %in% gwas_rnf144b_data$`ldbuddy_CHR:hg38POS`) |> 
  dplyr::select(`CHR:hg38POS`, chrom, hg38pos, p) |> 
  dplyr::rename(snp = `CHR:hg38POS`, pos = hg38pos) |> 
  mutate(ldbuddy_R2 = NA)

rnf144b_gwas_signal <- gwas_rnf144b_data |> 
  dplyr::select(-`ldbuddy_CHR:hg38POS`) |> 
  bind_rows(rnf144b_summary_stats) |> 
  mutate(chrom = paste0("chr", chrom)) |> 
  # Create group for LD group
  mutate(LDgrp = factor(cut(ldbuddy_R2, c(0, 0.2, 0.4, 0.6, 0.8, 1)),
                        levels = c("(0,0.2]", "(0.2,0.4]", "(0.4,0.6]",
                                   "(0.6,0.8]", "(0.8,1]", NA)))
rnf144b_gwas_signal$LDgrp <- addNA(rnf144b_gwas_signal$LDgrp)

# eQTL
rnf144b_signals <- bind_rows(fread(paste0("data/eqtl/qtl_nom/CTL_PEER_k20_genoPC_allSignals_nom1Mb_MAFs_", 
                                          rnf144b_egene[["gene_chr"]], ".csv"),
                                   data.table = FALSE) |> 
                               filter(gene_id == rnf144b_egene[["gene_id"]] & signal == rnf144b_egene[["signal"]]) |> 
                               mutate(Condition = "PBS"),
                             fread(paste0("data/eqtl/qtl_nom/FNF_PEER_k22_genoPC_allSignals_nom1Mb_MAFs_", 
                                          rnf144b_egene[["gene_chr"]], ".csv"),
                                   data.table = FALSE) |> 
                               filter(gene_id == rnf144b_egene[["gene_id"]] & signal == rnf144b_egene[["signal"]]) |> 
                               mutate(Condition = "FN-f")) |>
  mutate(rsID = ifelse(variantID == rnf144b_egene[["variantID"]], rnf144b_qtl_rsid, 
                       ifelse(variantID %in% rnf144b_gwas_lead_ids, rnf144b_gwas_rsid, NA))) |> 
  dplyr::rename(chrom = variant_chr,
                pos = variant_start,
                p = nom_pval,
                snp = rsID) |> 
  left_join(rnf144b_egene_ld, by = join_by(variantID == ld_variantID)) |> 
  mutate(LDgrp = cut(R2, c(0, 0.2, 0.4, 0.6, 0.8, 1)))
rnf144b_signals$LDgrp <- addNA(rnf144b_signals$LDgrp)

rnf144b_pbs_qval <- read_csv("data/eqtl/qtl_supp/CTL_PEER_k20_genoPC_perm1Mb_FDR.csv") |> 
  filter(gene_symbol == "RNF144B") |> 
  pull(qval)

rnf144b_fnf_qval <- read_csv("data/eqtl/qtl_supp/FNF_PEER_k22_genoPC_perm1Mb_FDR.csv") |> 
  filter(gene_symbol == "RNF144B") |> 
  pull(qval)
# Create figure with plotgardener -----------------------------------------

pdf("plots/colocs_Fig5/Fig5.pdf", width = 12, height = 7.5)
pageCreate(width = 12, height = 7.5, showGuides = FALSE)

## A- TGFA
plotText("A", x = 0.1, y = 0.1, just = c("left", "top"), fontfamily = "Helvetica",
         fontsize = 14, fontface = "bold")
plotGG(coloc_shared_eqtl_boxplot, x = 0.25, y = 0, width = 3.6, height = 2.85)
plotText(paste0("qval = ", signif(tgfa_pbs_qval, digits = 3)), 
         x = 1.98, y = 0.625, just = "right", fontfamily = "Helvetica",
         fontsize = 5)
plotText(paste0("qval = ", signif(tgfa_fnf_qval, digits = 3)), 
         x = 3.67, y = 0.625, just = "right", fontfamily = "Helvetica",
         fontsize = 5)


tgfa_region <- pgParams(assembly = "hg38",
                   chrom = tgfa_egene[["gene_chr"]],
                   chromstart = as.numeric(tgfa_egene[["variant_start"]]) - 500000,
                   chromend = as.numeric(tgfa_egene[["variant_start"]]) + 500000,
                   x = 0.5, width = 3.25)
qtl_tgfa_ylim <- ceiling(max(-1*log10(tgfa_signals$p))) + 1.5
gwas_tgfa_ylim <- ceiling(max(-1*log10(tgfa_gwas_signal$p))) + 1.5

gwas_tgfa <- plotManhattan(data = tgfa_gwas_signal |> 
                              arrange(desc(LDgrp)),
                            params = tgfa_region,
                            range = c(0, gwas_tgfa_ylim),
                            y = 2.75, height = 1.15, 
                            sigline = FALSE,
                            fill = colorby("LDgrp",
                                           palette = colorRampPalette(c("#262C74",
                                                                        "#98CDED",
                                                                        "#499A53",
                                                                        "#EEA741",
                                                                        "#DD3931",
                                                                        "grey"))),
                            snpHighlights = data.frame(snp = c(tgfa_qtl_rsid, tgfa_gwas_rsid),
                                                       pch = c(24, 23),
                                                       cex = c(0.75, 0.75),
                                                       col = c("black", "black")))

annoYaxis(plot = gwas_tgfa, at = seq(0, gwas_tgfa_ylim, 4), 
          axisLine = TRUE, fontsize = 7)
plotText(
  label = "-log10(p-value)", x = 0.15, y = 3.3, rot = 90,
  fontsize = 7, just = "center",
  default.units = "inches", fontfamily = "Helvetica"
)
plotText(label = tgfa_coloc$OA, x = 3.75, y = 2.75, just = "right", fontface = "bold",
         fontsize = 10, fontfamily = "Helvetica")
pbs_tgfa <- plotManhattan(data = tgfa_signals |> filter(Condition == "PBS"),
                           params = tgfa_region, 
                           range = c(0, qtl_tgfa_ylim),
                           snpHighlights = data.frame(snp = c(tgfa_gwas_rsid, 
                                                              tgfa_qtl_rsid),
                                                      pch = c(23, 24),
                                                      cex = c(0.75, 0.75),
                                                      col = c("black", "black")),
                             fill = colorby("LDgrp",
                                            palette = colorRampPalette(c("#262C74",
                                                                         "#98CDED",
                                                                         "#499A53",
                                                                         "#EEA741",
                                                                         "#DD3931",
                                                                         "grey"))),
                             y = 4, height = 1.15)


annoYaxis(plot = pbs_tgfa, at = seq(0, qtl_tgfa_ylim, 2), 
          axisLine = TRUE, fontsize = 7)
plotText(
  label = "-log10(p-value)", x = 0.15, y = 4.6, rot = 90,
  fontsize = 7, just = "center",
  default.units = "inches", fontfamily = "Helvetica"
)
plotText(label = "PBS", x = 3.75, y = 4.1, just = "right", fontface = "bold",
         fontsize = 10, fontfamily = "Helvetica")

fnf_tgfa <- plotManhattan(data = tgfa_signals |> filter(Condition == "FN-f"),
                           params = tgfa_region, 
                           range = c(0, qtl_tgfa_ylim),
                           snpHighlights = data.frame(snp = c(tgfa_gwas_rsid, 
                                                              tgfa_qtl_rsid),
                                                      pch = c(23, 24),
                                                      cex = c(0.75, 0.75),
                                                      col = c("black", "black")),
                           fill = colorby("LDgrp",
                                          palette = colorRampPalette(c("#262C74",
                                                                       "#98CDED",
                                                                       "#499A53",
                                                                       "#EEA741",
                                                                       "#DD3931",
                                                                       "grey"))),
                           y = 5.25, height = 1.15)

annoYaxis(plot = fnf_tgfa, at = seq(0, qtl_tgfa_ylim, 2), 
          axisLine = TRUE, fontsize = 7)
plotText(
  label = "-log10(p-value)", x = 0.15, y = 5.78, rot = 90,
  fontsize = 7, just = "center",
  default.units = "inches", fontfamily = "Helvetica"
)
plotText(label = "FN-f", x = 3.75, y = 5.4, just = "right", fontface = "bold",
         fontsize = 10, fontfamily = "Helvetica")

tgfa_genes <- plotGenes(params = tgfa_region, y = 6.5,
                   height = 0.5, geneHighlights = data.frame("gene" = "TGFA",
                                                             "color" = "#37a7db"))
annoGenomeLabel(plot = tgfa_genes, params = tgfa_region, y = 7, fontsize = 8)

plotText(label = paste0("PP4 = ", 
                        signif(tgfa_coloc |> pull(PP4_pbs), digits = 3)), 
         x = 3.75, y = 4.25, 
         just = "right", fontsize = 9, fontfamily = "Helvetica")

plotText(label = paste0("PP4 = ", 
                        signif(tgfa_coloc |> pull(PP4_fnf), digits = 3)), 
         x = 3.75, y = 5.55, 
         just = "right", fontsize = 9, fontfamily = "Helvetica")

grid.points(x = 0.95, y = 7.3, default.units = "native", pch = 23, size = unit(0.75, "char"))
plotText(label = paste0(tgfa_gwas_rsid, " (", tgfa_coloc$OA, " index)"),
         fontsize = 7, fontfamily = "Helvetica",
         just = "right", x = 2.1, y = 7.3)
grid.points(x = 2.25, y = 7.31, default.units = "native", pch = 24, size = unit(0.65, "char"))
plotText(label = paste0(tgfa_qtl_rsid, " (eQTL index)"),
         fontsize = 7, fontfamily = "Helvetica",
         just = "right", x = 3.4, y = 7.3)


plotLegend(legend = c("0.8 - 1.0",
                      "0.6 - 0.8",
                      "0.4 - 0.6",
                      "0.2 - 0.4",
                      "0.0 - 0.2"),
           fill = c("#DD3931", "#EEA741", "#499A53","#98CDED","#262C74"),
           x = 0.4, y = 2.9, width = 0.1, height = 0.4, border = FALSE, 
           fontsize = 6)
plotText(label = "r2", x = 1, y = 2.8, fontsize = 6, fontfamily = "Helvetica")

## B - PIK3R1
plotText("B", x = 4, y = 0.1, just = c("left", "top"), fontfamily = "Helvetica",
         fontsize = 14, fontface = "bold")
plotGG(coloc_pbs_eqtl_boxplot, x = 4.25, y = 0, width = 3.6, height = 2.85)

plotText(paste0("qval = ", signif(pik3r1_pbs_qval, digits = 3)), 
         x = 5.98, y = 0.625, just = "right", fontfamily = "Helvetica",
         fontsize = 5)
plotText(paste0("qval = ", signif(pik3r1_fnf_qval, digits = 3)), 
         x = 7.67, y = 0.625, just = "right", fontfamily = "Helvetica",
         fontsize = 5)
pik3r1_region <- pgParams(assembly = "hg38",
                         chrom = pik3r1_egene[["gene_chr"]],
                         chromstart = as.numeric(pik3r1_egene[["variant_start"]]) - 500000,
                         chromend = as.numeric(pik3r1_egene[["variant_start"]]) + 500000,
                         x = 4.5, width = 3.25)
qtl_pik3r1_ylim <- ceiling(max(-1*log10(pik3r1_signals$p))) + 1.5
gwas_pik3r1_ylim <- ceiling(max(-1*log10(pik3r1_gwas_signal$p))) + 1.5

gwas_pik3r1 <- plotManhattan(data = pik3r1_gwas_signal |> 
                              arrange(desc(LDgrp)),
                            params = pik3r1_region,
                            range = c(0, gwas_pik3r1_ylim),
                            y = 2.75, height = 1.15, 
                            sigline = FALSE,
                            fill = colorby("LDgrp",
                                           palette = colorRampPalette(c("#262C74",
                                                                        "#98CDED",
                                                                        "#499A53",
                                                                        "#EEA741",
                                                                        "#DD3931",
                                                                        "grey"))),
                            snpHighlights = data.frame(snp = c(pik3r1_qtl_rsid, pik3r1_gwas_rsid),
                                                       pch = c(23, 24),
                                                       cex = c(0.75, 0.75),
                                                       col = c("black", "black")))
annoYaxis(plot = gwas_pik3r1, at = seq(0, gwas_pik3r1_ylim, 2), 
          axisLine = TRUE, fontsize = 7)
plotText(
  label = "-log10(p-value)", x = 4.2, y = 3.3, rot = 90,
  fontsize = 7, just = "center",
  default.units = "inches", fontfamily = "Helvetica"
)
plotText(label = pik3r1_coloc$OA, x = 7.75, y = 2.75, just = "right", fontface = "bold",
         fontsize = 10, fontfamily = "Helvetica")

pbs_pik3r1 <- plotManhattan(data = pik3r1_signals |> filter(Condition == "PBS"),
                             params = pik3r1_region,
                             range = c(0, qtl_pik3r1_ylim),
                             y = 4, height = 1.15, 
                             sigline = FALSE,
                             fill = colorby("LDgrp",
                                            palette = colorRampPalette(c("#262C74",
                                                                         "#98CDED",
                                                                         "#499A53",
                                                                         "#EEA741",
                                                                         "#DD3931",
                                                                         "grey"))),
                             snpHighlights = data.frame(snp = c(pik3r1_gwas_rsid, pik3r1_qtl_rsid),
                                                        pch = c(24, 23),
                                                        cex = c(0.75, 0.75),
                                                        col = c("black", "black")))

annoYaxis(plot = pbs_pik3r1, at = seq(0, qtl_pik3r1_ylim, 1), 
          axisLine = TRUE, fontsize = 7)
plotText(
  label = "-log10(p-value)", x = 4.2, y = 4.6, rot = 90,
  fontsize = 7, just = "center",
  default.units = "inches", fontfamily = "Helvetica"
)
plotText(label = "PBS", x = 7.75, y = 4.1, just = "right", fontface = "bold",
         fontsize = 10, fontfamily = "Helvetica")

fnf_pik3r1 <- plotManhattan(data = pik3r1_signals |> filter(Condition == "FN-f"),
                            params = pik3r1_region,
                            range = c(0, qtl_pik3r1_ylim),
                            y = 5.25, height = 1.15, 
                            sigline = FALSE,
                            fill = colorby("LDgrp",
                                           palette = colorRampPalette(c("#262C74",
                                                                        "#98CDED",
                                                                        "#499A53",
                                                                        "#EEA741",
                                                                        "#DD3931",
                                                                        "grey"))),
                            snpHighlights = data.frame(snp = c(pik3r1_gwas_rsid, pik3r1_qtl_rsid),
                                                       pch = c(24, 23),
                                                       cex = c(0.75, 0.75),
                                                       col = c("black", "black")))

annoYaxis(plot = fnf_pik3r1, at = seq(0, qtl_pik3r1_ylim, 1), 
          axisLine = TRUE, fontsize = 7)
plotText(
  label = "-log10(p-value)", x = 4.2, y = 5.78, rot = 90,
  fontsize = 7, just = "center",
  default.units = "inches", fontfamily = "Helvetica"
)
plotText(label = "FN-f", x = 7.75, y = 5.4, just = "right", fontface = "bold",
         fontsize = 10, fontfamily = "Helvetica")

pik3r1_genes <- plotGenes(params = pik3r1_region, y = 6.5,
                         height = 0.5, geneHighlights = data.frame("gene" = "PIK3R1",
                                                                   "color" = "#37a7db"))
annoGenomeLabel(plot = pik3r1_genes, params = pik3r1_region, y = 7, fontsize = 8)

plotText(label = paste0("PP4 = ", 
                        signif(pik3r1_coloc |> pull(PP4_pbs), digits = 3)), 
         x = 7.75, y = 4.25, 
         just = "right", fontsize = 9, fontfamily = "Helvetica")

plotText(label = paste0("PP4 = ", 
                        signif(pik3r1_coloc |> pull(PP4_fnf), digits = 3)), 
         x = 7.75, y = 5.55, 
         just = "right", fontsize = 9, fontfamily = "Helvetica")

grid.points(x = 5, y = 7.31, default.units = "native", pch = 24, 
            size = unit(0.65, "char"))
plotText(label = paste0(pik3r1_gwas_rsid, " (", pik3r1_coloc$OA, " index)"),
         fontsize = 7, fontfamily = "Helvetica",
         just = "left", x = 5.1, y = 7.3)

grid.points(x = 6.3, y = 7.3, default.units = "native", pch = 23, 
            size = unit(0.75, "char"))
plotText(label = paste0(pik3r1_qtl_rsid, " (eQTL index)"),
         fontsize = 7, fontfamily = "Helvetica",
         just = "left", x = 6.4, y = 7.3)

plotLegend(legend = c("0.8 - 1.0",
                      "0.6 - 0.8",
                      "0.4 - 0.6",
                      "0.2 - 0.4",
                      "0.0 - 0.2"),
           fill = c("#DD3931", "#EEA741", "#499A53","#98CDED","#262C74"),
           x = 4.4, y = 2.9, width = 0.1, height = 0.4, border = FALSE, 
           fontsize = 6)
plotText(label = "r2", x = 5, y = 2.8, fontsize = 6, fontfamily = "Helvetica")

## C - RNF144B
plotText("C", x = 8, y = 0.1, just = c("left", "top"), fontfamily = "Helvetica",
         fontsize = 14, fontface = "bold")
plotGG(coloc_fnf_eqtl_boxplot, x = 8.15, y = 0, width = 3.6, height = 2.85)

plotText(paste0("qval = ", signif(rnf144b_pbs_qval, digits = 3)), 
         x = 9.88, y = 0.625, just = "right", fontfamily = "Helvetica",
         fontsize = 5)
plotText(paste0("qval = ", signif(rnf144b_fnf_qval, digits = 3)), 
         x = 11.57, y = 0.625, just = "right", fontfamily = "Helvetica",
         fontsize = 5)

rnf144b_region <- pgParams(assembly = "hg38",
                          chrom = rnf144b_egene[["gene_chr"]],
                          chromstart = as.numeric(rnf144b_egene[["variant_start"]]) - 500000,
                          chromend = as.numeric(rnf144b_egene[["variant_start"]]) + 500000,
                          x = 8.5, width = 3.25)
qtl_rnf144b_ylim <- ceiling(max(-1*log10(rnf144b_signals$p))) + 1.5
gwas_rnf144b_ylim <- ceiling(max(-1*log10(rnf144b_gwas_signal$p))) + 1.5

gwas_rnf144b <- plotManhattan(data = rnf144b_gwas_signal |> 
                               arrange(desc(LDgrp)),
                             params = rnf144b_region,
                             range = c(0, gwas_rnf144b_ylim),
                             y = 2.75, height = 1.15, 
                             sigline = FALSE,
                             fill = colorby("LDgrp",
                                            palette = colorRampPalette(c("#262C74",
                                                                         "#98CDED",
                                                                         "#499A53",
                                                                         "#EEA741",
                                                                         "#DD3931",
                                                                         "grey"))),
                             snpHighlights = data.frame(snp = c(rnf144b_qtl_rsid, rnf144b_gwas_rsid),
                                                        pch = c(23, 24),
                                                        cex = c(0.75, 0.75),
                                                        col = c("black", "black")))
annoYaxis(plot = gwas_rnf144b, at = seq(0, gwas_rnf144b_ylim, 2), 
          axisLine = TRUE, fontsize = 7)
plotText(
  label = "-log10(p-value)", x = 8.2, y = 3.3, rot = 90,
  fontsize = 7, just = "center",
  default.units = "inches", fontfamily = "Helvetica"
)
plotText(label = rnf144b_coloc$OA, x = 11.75, y = 2.75, just = "right", fontface = "bold",
         fontsize = 10, fontfamily = "Helvetica")

pbs_rnf144b <- plotManhattan(data = rnf144b_signals |> filter(Condition == "PBS"),
                            params = rnf144b_region,
                            range = c(0, qtl_rnf144b_ylim),
                            y = 4, height = 1.15, 
                            sigline = FALSE,
                            fill = colorby("LDgrp",
                                           palette = colorRampPalette(c("#262C74",
                                                                        "#98CDED",
                                                                        "#499A53",
                                                                        "#EEA741",
                                                                        "#DD3931",
                                                                        "grey"))),
                            snpHighlights = data.frame(snp = c(rnf144b_gwas_rsid, rnf144b_qtl_rsid),
                                                       pch = c(24, 23),
                                                       cex = c(0.75, 0.75),
                                                       col = c("black", "black")))

annoYaxis(plot = pbs_rnf144b, at = seq(0, qtl_rnf144b_ylim, 1), 
          axisLine = TRUE, fontsize = 7)
plotText(
  label = "-log10(p-value)", x = 8.2, y = 4.6, rot = 90,
  fontsize = 7, just = "center",
  default.units = "inches", fontfamily = "Helvetica"
)
plotText(label = "PBS", x = 11.75, y = 4.1, just = "right", fontface = "bold",
         fontsize = 10, fontfamily = "Helvetica")

fnf_rnf144b <- plotManhattan(data = rnf144b_signals |> filter(Condition == "FN-f"),
                            params = rnf144b_region,
                            range = c(0, qtl_rnf144b_ylim),
                            y = 5.25, height = 1.15, 
                            sigline = FALSE,
                            fill = colorby("LDgrp",
                                           palette = colorRampPalette(c("#262C74",
                                                                        "#98CDED",
                                                                        "#499A53",
                                                                        "#EEA741",
                                                                        "#DD3931",
                                                                        "grey"))),
                            snpHighlights = data.frame(snp = c(rnf144b_gwas_rsid, rnf144b_qtl_rsid),
                                                       pch = c(24, 23),
                                                       cex = c(0.75, 0.75),
                                                       col = c("black", "black")))

annoYaxis(plot = fnf_rnf144b, at = seq(0, qtl_rnf144b_ylim, 1), 
          axisLine = TRUE, fontsize = 7)
plotText(
  label = "-log10(p-value)", x = 8.2, y = 5.78, rot = 90,
  fontsize = 7, just = "center",
  default.units = "inches", fontfamily = "Helvetica"
)
plotText(label = "FN-f", x = 11.75, y = 5.4, just = "right", fontface = "bold",
         fontsize = 10, fontfamily = "Helvetica")

rnf144b_genes <- plotGenes(params = rnf144b_region, y = 6.5,
                          height = 0.5, geneHighlights = data.frame("gene" = "RNF144B",
                                                                    "color" = "#37a7db"))
plotText(label = paste0("PP4 = ", 
                        signif(rnf144b_coloc |> pull(PP4_pbs), digits = 3)), 
         x = 11.75, y = 4.25, 
         just = "right", fontsize = 9, fontfamily = "Helvetica")
annoGenomeLabel(plot = rnf144b_genes, params = rnf144b_region, y = 7, fontsize = 8)

plotText(label = paste0("PP4 = ", 
                        signif(rnf144b_coloc |> pull(PP4_fnf), digits = 3)), 
         x = 11.75, y = 5.55, 
         just = "right", fontsize = 9, fontfamily = "Helvetica")

grid.points(x = 8.9, y = 7.31, default.units = "native", pch = 24, 
            size = unit(0.65, "char"))
plotText(label = paste0(rnf144b_gwas_rsid, " (", rnf144b_coloc$OA, " index)"),
         fontsize = 7, fontfamily = "Helvetica",
         just = "left", x = 9, y = 7.3)

grid.points(x = 10.35, y = 7.3, default.units = "native", pch = 23, 
            size = unit(0.75, "char"))
plotText(label = paste0(rnf144b_qtl_rsid, " (eQTL index)"),
         fontsize = 7, fontfamily = "Helvetica",
         just = "left", x = 10.45, y = 7.3)

plotLegend(legend = c("0.8 - 1.0",
                      "0.6 - 0.8",
                      "0.4 - 0.6",
                      "0.2 - 0.4",
                      "0.0 - 0.2"),
           fill = c("#DD3931", "#EEA741", "#499A53","#98CDED","#262C74"),
           x = 8.4, y = 2.9, width = 0.1, height = 0.4, border = FALSE, 
           fontsize = 6)
plotText(label = "r2", x = 9, y = 2.8, fontsize = 6, fontfamily = "Helvetica")
dev.off()
