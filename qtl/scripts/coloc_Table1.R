library(tidyverse)
library(googledrive)
library(googlesheets4)

# Wrangle gene/qtl data ---------------------------------------------------

# condition-specific reGenes
reGenes <- full_join(read_csv("data/reqtl/CTL_sig01_reQTLs_PEER_k20_genoPC.csv", 
                              col_select = c("gene_id")) |> 
                       mutate(pbs_reGene = "PBS"),
                     read_csv("data/reqtl/FNF_sig01_reQTLs_PEER_k22_genoPC.csv", 
                              col_select = c("gene_id")) |> 
                       mutate(fnf_reGene = "FN-f"), by = "gene_id") |> 
  rowwise() |> 
  mutate(reGene = paste(unique(na.omit(c(pbs_reGene, fnf_reGene))), collapse = "/")) |> 
  ungroup() |> 
  mutate(reGene = ifelse(reGene == "PBS/FN-f", "Both", reGene)) |> 
  dplyr::select(-pbs_reGene, -fnf_reGene)

# Overlap with OA, multiple studies
raak_de_genes <- read_csv("../DE/data/RAAK/RAAK_genes.csv", col_select = c("HGNC", 
                                                                           "RAAK_LFC")) |> 
  mutate(raak_dir = ifelse(RAAK_LFC < 0, "down", "up")) |> 
  dplyr::rename(symbol = HGNC) |> 
  dplyr::select(-RAAK_LFC)

fisch2018_de_genes <- read_csv("../DE/data/GSE114007/1-s2.0-S1063458418313876-mmc1.csv",
                               col_names = c("symbol", "log2FoldChange", "padj"),
                               col_select = c(1,2, 3), n_max = 12475) |> 
  filter(padj < 0.05) |> 
  mutate(fisch_dir = ifelse(log2FoldChange > 0, "up", "down")) |> 
  dplyr::select(-log2FoldChange, -padj)


fu2021_de_genes <- read_csv("../DE/data/GSE168505/GSE168505_deseq_res.csv", 
                            col_select = c("symbol", "log2FoldChange", "padj")) |> 
  filter(padj < 0.05) |> 
  mutate(fu_dir = ifelse(log2FoldChange < 0, "down", "up")) |> 
  dplyr::select(-log2FoldChange, -padj)

oa_de_genes <- full_join(raak_de_genes, fisch2018_de_genes, by = "symbol") |> 
  full_join(fu2021_de_genes, by = "symbol") |> 
  rowwise() |> 
  mutate(oa_dir = paste(unique(na.omit(c(raak_dir, fisch_dir, fu_dir))), collapse = "/")) |> 
  ungroup() |> 
  dplyr::select(-raak_dir, -fisch_dir, -fu_dir)

# Overlap with FN-f DE
fnf_de_genes <- read_csv("../DE/data/condition_de/sig_deGenes_pval01_l2fc1.csv",
                      col_select = c("gene_id", "log2FoldChange")) |> 
  mutate(fnf_dir = ifelse(log2FoldChange < 0, "down", "up")) |> 
  dplyr::select(-log2FoldChange)


# Overlap with age-related expression
age_de_genes <- full_join(read_csv("../DE/data/age_de/ctl_age_pval05clusters.csv",
                                   col_select = c("gene_id", "cluster")) |> 
                            distinct() |>
                            mutate(cluster = ifelse(cluster == "+", "up", "down")) |> 
                            dplyr::rename(pbs_age_cluster = cluster),
                          read_csv("../DE/data/age_de/fnf_age_pval05clusters.csv",
                                   col_select = c("gene_id", "cluster")) |> 
                            distinct() |> 
                            mutate(cluster = ifelse(cluster == "+", "up", "down")) |> 
                            dplyr::rename(fnf_age_cluster = cluster), by = "gene_id") |> 
  rowwise() |> 
  mutate(age_dir = paste(unique(na.omit(c(pbs_age_cluster, fnf_age_cluster))), collapse = "/")) |> 
  ungroup() |> 
  dplyr::select(-pbs_age_cluster, -fnf_age_cluster)


# Overlap with sex DE
sex_de_genes <- full_join(read_csv("../DE/data/sex_de/ctl_sexDE_pval01.csv", 
                      col_select = c("gene_id", "log2FoldChange")) |> 
                        mutate(pbs_sex_dir = ifelse(log2FoldChange < 0, "female", "male")) |> 
                        dplyr::select(-log2FoldChange),
read_csv("../DE/data/sex_de/fnf_sexDE_pval01.csv", 
         col_select = c("gene_id", "log2FoldChange")) |> 
  mutate(fnf_sex_dir = ifelse(log2FoldChange < 0, "female", "male")) |> 
  dplyr::select(-log2FoldChange), by = "gene_id") |> 
  rowwise() |> 
  mutate(sex_dir = paste(unique(na.omit(c(pbs_sex_dir, fnf_sex_dir))), collapse = "/")) |> 
  ungroup() |> 
  dplyr::select(-pbs_sex_dir, -fnf_sex_dir)
# Join all summaries into table -------------------------------------------

# GWAS/eQTL formatting
## All columns for supplement
coloc_gwas_qtl_gene_table <- read_csv("data/colocalization/processed/PBS_FNF_allSignals_EUR_GWAS_250Kb_colocs.csv",
         col_select = c("eGene_id", "eGene_name", 
                        "OA", "GWAS_lead", "GWAS_lead_variantID", "GWAS_OR", "GWAS_risk_allele", "GWAS_risk_allele_frequency",
                        "Condition", "qtl_rsid", "qtl_variantID", "PP4_pbs", "PP4_fnf", "riskAllele_eqtl_beta")) |> 
  separate_wider_delim(cols = "GWAS_lead_variantID", delim = ":", names = c("GWAS_chr", "GWAS_pos", "GWAS_A1", "GWAS_A2")) |> 
  separate_wider_delim(cols = "qtl_variantID", delim = ":", names = c("qtl_chr", "qtl_pos", NA, NA)) |> 
  mutate(OA = gsub("(.*)OA$", "\\1 OA", `OA`),
         `GWAS lead SNP position` = paste0(GWAS_chr, ":", GWAS_pos),
         `protective allele` = ifelse(GWAS_risk_allele == GWAS_A1, GWAS_A2, GWAS_A1),
         `eQTL lead SNP position` = paste0(qtl_chr, ":", qtl_pos),
         `risk allele association with expression` = ifelse(riskAllele_eqtl_beta < 0, "lower", "higher"),
        `novel colocalization` = ifelse(eGene_name %in% c("NPC1", "FAM53A", "SMAD3",
                                                   "SLC44A2", "ALDH1A2"), "no", "yes")) |> 
  dplyr::rename(Gene = eGene_name, `GWAS top OA phenotype` = OA,
                `GWAS OR` = GWAS_OR,
                `GWAS lead SNP` = GWAS_lead,
                `risk allele` = GWAS_risk_allele,
                `risk allele frequency` = GWAS_risk_allele_frequency,
                `eQTL condition` = Condition,
                `eQTL lead SNP` = qtl_rsid,
                `PBS coloc PP4` = PP4_pbs,
                `FN-f coloc PP4` = PP4_fnf) |> 
  dplyr::select(-GWAS_chr, -GWAS_pos, -GWAS_A1, -GWAS_A2, -riskAllele_eqtl_beta, -qtl_chr, -qtl_pos) |> 
  group_by(Gene) |> 
  mutate(`eQTL condition` = paste(`eQTL condition`, collapse = "/"),
         `eQTL lead SNP` = paste(`eQTL lead SNP`, collapse = "/"),
          `eQTL lead SNP position` = paste(`eQTL lead SNP position`, collapse = "/"),
         `PBS coloc PP4` = paste(`PBS coloc PP4`, collapse = "/"),
         `FN-f coloc PP4` = paste(`FN-f coloc PP4`, collapse = "/")) |> 
  ungroup() |> 
  distinct() |> 
  mutate(`eQTL condition` = ifelse(`eQTL condition` == "PBS/FN-f", "Both", `eQTL condition`)) |> 
  ## Joining with other gene/qtl data
  # Condition specific eQTL
  left_join(reGenes, by = join_by("eGene_id" == "gene_id")) |> 
  mutate(reGene = replace_na(reGene, "-")) |> 
  dplyr::rename(`eQTL genetic and condition interaction` = reGene) |> 
  # OA
  left_join(oa_de_genes, by = join_by("Gene" == "symbol")) |> 
  mutate(oa_dir = replace_na(oa_dir, "-")) |> 
  dplyr::rename(`change in OA` = oa_dir) |> 
  # FN-f
  left_join(fnf_de_genes, by = join_by("eGene_id" == "gene_id")) |> 
  mutate(fnf_dir = replace_na(fnf_dir, "-")) |> 
  dplyr::rename(`change with FN-f` = fnf_dir) |> 
  # Age
  left_join(age_de_genes, by = join_by("eGene_id" == "gene_id")) |> 
  mutate(age_dir = replace_na(age_dir, "-")) |>
  dplyr::rename(`change with age` = age_dir) |> 
  # Sex
  left_join(sex_de_genes, by = join_by("eGene_id" == "gene_id")) |> 
  mutate(sex_dir = replace_na(sex_dir, "-")) |>
  dplyr::rename(`M vs F` = sex_dir) |> 
  ## Clean up columns and reorder
  dplyr::select(-eGene_id) |> 
  relocate(`GWAS lead SNP position`, .after = `GWAS lead SNP`) |> 
  relocate(`protective allele`, .after = `risk allele`) |> 
  relocate(`risk allele association with expression`, .after = `risk allele frequency`) |> 
  relocate(`eQTL genetic and condition interaction`, .after = `eQTL condition`) |> 
  relocate(`eQTL lead SNP position`, .after = `eQTL lead SNP`) |> 
  arrange(Gene)

# Write to csv
write_csv(coloc_gwas_qtl_gene_table, file = "tables/SupTable13.csv")


# Write to google drive

ss <- gs4_create(name = "SupTable13")

write_sheet(coloc_gwas_qtl_gene_table,
            ss, sheet = "Sheet1")

drive_mv(file = "SupTable13", path = as_dribble("chQTL/chQTL paper/Figures and Tables"))

## Trim columns for main
coloc_gwas_qtl_gene_table_main <- coloc_gwas_qtl_gene_table |> 
  rowwise() |> 
  mutate(PP4 = case_when(`eQTL condition` == "PBS" ~ as.character(signif(as.numeric(`PBS coloc PP4`), digits = 3)),
                         `eQTL condition` == "FN-f" ~ as.character(signif(as.numeric(`FN-f coloc PP4`), digits = 3)),
                         `eQTL condition` == "Both" ~ paste0(signif(as.numeric(unlist(str_split(`PBS coloc PP4`, "/"))[1]), 
                                                                    digits = 3),
                                                             "/",
                                                             signif(as.numeric(unlist(str_split(`FN-f coloc PP4`, "/"))[2]), 
                                                                    digits = 3)))) |> 
  ungroup() |> 
  dplyr::select(-`GWAS lead SNP position`, -`GWAS OR`, 
                -`eQTL lead SNP position`, -`protective allele`,
                -`eQTL lead SNP`, -`PBS coloc PP4`, -`FN-f coloc PP4`) |> 
  relocate(PP4, .after = "eQTL condition") |> 
  relocate(`novel colocalization`, .after = "PP4")


# Write to csv
write_csv(coloc_gwas_qtl_gene_table_main, file = "tables/Table1.csv")


# Write to google drive

ss <- gs4_create(name = "Table1")

write_sheet(coloc_gwas_qtl_gene_table_main,
            ss, sheet = "Sheet1")

drive_mv(file = "Table1", path = as_dribble("chQTL/chQTL paper/Figures and Tables"))
