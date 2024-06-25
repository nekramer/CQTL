library(tidyverse)
library(plotgardener)
library(vcfR)
library(patchwork)
source("../plotting_utils.R")
source("../utils.R")



# Boxplots of GWAS variants shown for PBS and FN-f ------------------------

## TGFA

tgfa_coloc <- read_csv("data/colocalization/processed/PBS_FNF_allSignals_EUR_GWAS_250Kb_colocs.csv") |> 
  filter(eGene_name == "TGFA") |> 
  dplyr::select(GWAS_lead_variantID, eGene_id, signal, eGene_name, GWAS_lead, GWAS_risk_allele) |> 
  rename(variantID = GWAS_lead_variantID,
         gene_id = eGene_id,
         gene_symbol = eGene_name,
         rsID = GWAS_lead,
         risk_allele = GWAS_risk_allele)

# Get minor allele for genotype ordering 
tgfa_coloc$MA <- get_minor_allele(varID = tgfa_coloc$variantID, 
                                        data = tgfa_coloc)

# Get PBS and FNF betas and nominal p-vals
tgfa_egene_gwas <- tgfa_coloc |> 
  bind_cols(get_pvals_betas(varID = tgfa_coloc$variantID, data = tgfa_coloc))


vcf_chr2 <- vcfR2tidy(read.vcfR(paste0("data/vcf/CQTL_COA_01_GDA8_COA2_01_COA3_01_GDA8_COA4_COA5_COA6_COA7_qtl_chr2.vcf.gz"), 
                           verbose = FALSE))
variants_chr2 <- vcf_chr2$fix |> 
  dplyr::select(POS, ID, REF, ALT)
tgfa_geno_data <- vcf_chr2$gt |> 
  left_join(variants_chr2, by = "POS", relationship = "many-to-many") |> 
  filter(ID == tgfa_coloc$variantID) |> 
  dplyr::select(Indiv, gt_GT_alleles, ID) |> 
  dplyr::rename(Donor = Indiv,
                variantID = ID)

# Normalized expression
CTL_normQuant_tgfa <- read_delim("data/rna/CTL_CPMadjTMM_invNorm.bed.gz") |> 
  filter(gene_id == tgfa_coloc$gene_id) |> 
  dplyr::select(-`#chr`, -start, -end, -length, -strand) |> 
  pivot_longer(cols = starts_with("AM"), names_to = "Donor", values_to = "expression") |> 
  mutate(Condition = "PBS")

FNF_normQuant_tgfa <- read_delim("data/rna/FNF_CPMadjTMM_invNorm.bed.gz") |> 
  filter(gene_id == tgfa_coloc$gene_id)  |> 
  dplyr::select(-`#chr`, -start, -end, -length, -strand) |> 
  pivot_longer(cols = starts_with("AM"), names_to = "Donor", values_to = "expression") |> 
  mutate(Condition = "FN-f")

ALL_normQuant_tgfa <- bind_rows(CTL_normQuant_tgfa, FNF_normQuant_tgfa)

# Join all data elements together
tgfa_boxplot_data_gwas <- left_join(tgfa_egene_gwas, tgfa_geno_data, by = "variantID") |> 
  left_join(ALL_normQuant_tgfa, by = c("gene_id", "Donor"), relationship = "many-to-many")


# Get genotype order 
geno_order_tgfa <- determine_geno_order(varID = tgfa_coloc$variantID, 
                                        data = tgfa_boxplot_data_gwas)

tgfa_boxplot_data_gwas <- tgfa_boxplot_data_gwas |> 
  left_join(geno_order_tgfa, by = c("variantID", "gt_GT_alleles"))
tgfa_boxplot_data_gwas$Condition <- factor(tgfa_boxplot_data_gwas$Condition, 
                                            levels = c("PBS", "FN-f"))

tgfa_boxplot_gwas <- create_eqtl_boxplot_horizontal(tgfa_boxplot_data_gwas,
                                                     stat_loc = "top",
                                                     highlightAllele = "risk",
                                                     condition_labs = "bottom")

## PIK3R1

pik3r1_coloc <- read_csv("data/colocalization/processed/PBS_FNF_allSignals_EUR_GWAS_250Kb_colocs.csv") |> 
  filter(eGene_name == "PIK3R1") |> 
  dplyr::select(GWAS_lead_variantID, eGene_id, signal, eGene_name, GWAS_lead, GWAS_risk_allele) |> 
  rename(variantID = GWAS_lead_variantID,
         gene_id = eGene_id,
         gene_symbol = eGene_name,
         rsID = GWAS_lead,
         risk_allele = GWAS_risk_allele)

# Get minor allele for correct genotype ordering 
pik3r1_coloc$MA <- get_minor_allele(varID = pik3r1_coloc$variantID, 
                                  data = pik3r1_coloc)

# Get PBS and FNF betas and nominal p-vals
pik3r1_egene_gwas <- pik3r1_coloc |> 
  bind_cols(get_pvals_betas(varID = pik3r1_coloc$variantID, data = pik3r1_coloc))

vcf_chr5 <- vcfR2tidy(read.vcfR(paste0("data/vcf/CQTL_COA_01_GDA8_COA2_01_COA3_01_GDA8_COA4_COA5_COA6_COA7_qtl_chr5.vcf.gz"), 
                                verbose = FALSE))

variants_chr5 <- vcf_chr5$fix |> 
  dplyr::select(POS, ID, REF, ALT)
pik3r1_geno_data <- vcf_chr5$gt |> 
  left_join(variants_chr5, by = "POS", relationship = "many-to-many") |> 
  filter(ID == pik3r1_coloc$variantID) |> 
  dplyr::select(Indiv, gt_GT_alleles, ID) |> 
  dplyr::rename(Donor = Indiv,
                variantID = ID)

# Normalized expression
CTL_normQuant_pik3r1 <- read_delim("data/rna/CTL_CPMadjTMM_invNorm.bed.gz") |> 
  filter(gene_id == pik3r1_coloc$gene_id) |> 
  dplyr::select(-`#chr`, -start, -end, -length, -strand) |> 
  pivot_longer(cols = starts_with("AM"), names_to = "Donor", values_to = "expression") |> 
  mutate(Condition = "PBS")

FNF_normQuant_pik3r1 <- read_delim("data/rna/FNF_CPMadjTMM_invNorm.bed.gz") |> 
  filter(gene_id == pik3r1_coloc$gene_id)  |> 
  dplyr::select(-`#chr`, -start, -end, -length, -strand) |> 
  pivot_longer(cols = starts_with("AM"), names_to = "Donor", values_to = "expression") |> 
  mutate(Condition = "FN-f")

ALL_normQuant_pik3r1 <- bind_rows(CTL_normQuant_pik3r1, FNF_normQuant_pik3r1)

# Join all data elements together
pik3r1_boxplot_data_gwas <- left_join(pik3r1_egene_gwas, pik3r1_geno_data, by = "variantID") |> 
  left_join(ALL_normQuant_pik3r1, by = c("gene_id", "Donor"), relationship = "many-to-many")


# Get genotype order 
geno_order_pik3r1 <- determine_geno_order(varID = pik3r1_coloc$variantID, 
                                        data = pik3r1_boxplot_data_gwas)

pik3r1_boxplot_data_gwas <- pik3r1_boxplot_data_gwas |> 
  left_join(geno_order_pik3r1, by = c("variantID", "gt_GT_alleles"))
pik3r1_boxplot_data_gwas$Condition <- factor(pik3r1_boxplot_data_gwas$Condition, 
                                           levels = c("PBS", "FN-f"))

pik3r1_boxplot_gwas <- create_eqtl_boxplot_horizontal(pik3r1_boxplot_data_gwas,
                                                    stat_loc = "top",
                                                    highlightAllele = "risk",
                                                    condition_labs = "bottom")

## RNF144B

rnf144b_coloc <- read_csv("data/colocalization/processed/PBS_FNF_allSignals_EUR_GWAS_250Kb_colocs.csv") |> 
  filter(eGene_name == "RNF144B") |> 
  dplyr::select(GWAS_lead_variantID, eGene_id, signal, eGene_name, GWAS_lead, GWAS_risk_allele) |> 
  rename(variantID = GWAS_lead_variantID,
         gene_id = eGene_id,
         gene_symbol = eGene_name,
         rsID = GWAS_lead,
         risk_allele = GWAS_risk_allele) |> 
  mutate(variantID = paste0(unlist(str_split(variantID, ":"))[1], ":",
                            unlist(str_split(variantID, ":"))[2], ":",
                            unlist(str_split(variantID, ":"))[4], ":",
                            unlist(str_split(variantID, ":"))[3]))

# Get minor allele for correct genotype ordering 
rnf144b_coloc$MA <- get_minor_allele(varID = rnf144b_coloc$variantID, 
                                    data = rnf144b_coloc)

# Get PBS and FNF betas and nominal p-vals
rnf144b_egene_gwas <- rnf144b_coloc |> 
  bind_cols(get_pvals_betas(varID = rnf144b_coloc$variantID, data = rnf144b_coloc))

vcf_chr6 <- vcfR2tidy(read.vcfR(paste0("data/vcf/CQTL_COA_01_GDA8_COA2_01_COA3_01_GDA8_COA4_COA5_COA6_COA7_qtl_chr6.vcf.gz"), 
                                verbose = FALSE))

variants_chr6 <- vcf_chr6$fix |> 
  dplyr::select(POS, ID, REF, ALT)
rnf144b_geno_data <- vcf_chr6$gt |> 
  left_join(variants_chr6, by = "POS", relationship = "many-to-many") |> 
  filter(ID == rnf144b_coloc$variantID) |> 
  dplyr::select(Indiv, gt_GT_alleles, ID) |> 
  dplyr::rename(Donor = Indiv,
                variantID = ID)

# Normalized expression
CTL_normQuant_rnf144b <- read_delim("data/rna/CTL_CPMadjTMM_invNorm.bed.gz") |> 
  filter(gene_id == rnf144b_coloc$gene_id) |> 
  dplyr::select(-`#chr`, -start, -end, -length, -strand) |> 
  pivot_longer(cols = starts_with("AM"), names_to = "Donor", values_to = "expression") |> 
  mutate(Condition = "PBS")

FNF_normQuant_rnf144b <- read_delim("data/rna/FNF_CPMadjTMM_invNorm.bed.gz") |> 
  filter(gene_id == rnf144b_coloc$gene_id)  |> 
  dplyr::select(-`#chr`, -start, -end, -length, -strand) |> 
  pivot_longer(cols = starts_with("AM"), names_to = "Donor", values_to = "expression") |> 
  mutate(Condition = "FN-f")

ALL_normQuant_rnf144b <- bind_rows(CTL_normQuant_rnf144b, FNF_normQuant_rnf144b)

# Join all data elements together
rnf144b_boxplot_data_gwas <- left_join(rnf144b_egene_gwas, rnf144b_geno_data, by = "variantID") |> 
  left_join(ALL_normQuant_rnf144b, by = c("gene_id", "Donor"), relationship = "many-to-many")

# Get genotype order 
geno_order_rnf144b<- determine_geno_order(varID = rnf144b_coloc$variantID, 
                                          data = rnf144b_boxplot_data_gwas)

rnf144b_boxplot_data_gwas <- rnf144b_boxplot_data_gwas |> 
  left_join(geno_order_rnf144b, by = c("variantID", "gt_GT_alleles"))
rnf144b_boxplot_data_gwas$Condition <- factor(rnf144b_boxplot_data_gwas$Condition, 
                                             levels = c("PBS", "FN-f"))

rnf144b_boxplot_gwas <- create_eqtl_boxplot_horizontal(rnf144b_boxplot_data_gwas,
                                                      stat_loc = "top",
                                                      highlightAllele = "risk",
                                                      condition_labs = "bottom")



coloc_gwas_boxplots <- tgfa_boxplot_gwas + pik3r1_boxplot_gwas + rnf144b_boxplot_gwas +
  plot_annotation(theme = theme(panel.background = 
                                  element_rect(fill = "transparent", 
                                               color = "transparent"),
                                plot.background = 
                                  element_rect(fill = "transparent", 
                                               color = "transparent")))
save(coloc_gwas_boxplots, file = "plots/colocs_Fig5_supp/coloc_gwas_boxplots.rda")


pdf("plots/colocs_Fig5_supp/SupFig8.pdf", width = 6.6, height = 2.1)
pageCreate(width = 6.6, height = 2.1, showGuides = FALSE)
plotText("A", x = 0.1, y = 0.1, just = c("left", "top"), fontfamily = "Helvetica",
         fontsize = 14, fontface = "bold")
plotText("B", x = 2.2, y = 0.1, just = c("left", "top"), fontfamily = "Helvetica",
         fontsize = 14, fontface = "bold")
plotText("C", x = 4.4, y = 0.1, just = c("left", "top"), fontfamily = "Helvetica",
         fontsize = 14, fontface = "bold")
plotGG(coloc_gwas_boxplots, x = 0.1, y = 0, width = 6.5, height = 2.3)
dev.off()
