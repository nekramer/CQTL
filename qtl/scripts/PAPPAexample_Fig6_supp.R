library(tidyverse)
library(patchwork)
library(plotgardener)
source("../utils.R")
source("../plotting_utils.R")


# eQTL boxplots of PAPPA signal lead in PBS and FN-f ----------------------

## eGene information from both conditions
pappa_egene_pbs <- read_csv("data/eqtl/CTL_PEER_k20_genoPC_perm1Mb_sig_rsID_signalRanges.csv") |> 
  filter(gene_symbol == "PAPPA") |> 
  mutate(grouping = "PBS")

pappa_egene_fnf <- read_csv("data/eqtl/FNF_PEER_k22_genoPC_perm1Mb_sig_rsID_signalRanges.csv") |> 
  filter(gene_symbol == "PAPPA") |> 
  mutate(grouping = "FN-f")

# Join
pappa_egenes_qtls <- bind_rows(pappa_egene_pbs, pappa_egene_fnf)

# Get minor allele for each lead
pappa_egenes_qtls$MA <- unlist(lapply(pappa_egenes_qtls$variantID, get_minor_allele, pappa_egenes_qtls))

# Get nominal p-values and betas for each lead in both conditions
pappa_egenes_qtls <- pappa_egenes_qtls |> 
  bind_cols(lapply(pappa_egenes_qtls$variantID, get_pvals_betas, pappa_egenes_qtls) |> 
              bind_rows())

## VCF for genotypes
vcf <- vcfR2tidy(read.vcfR("data/vcf/CTLk20_FNFk22_genoPC_leadVars.vcf.gz"), 
                           verbose = FALSE)
variants <- vcf$fix |> 
  dplyr::select(POS, ID, REF, ALT)
geno_data <- vcf$gt |> 
  left_join(variants, by = "POS", relationship = "many-to-many") |> 
  filter(ID %in% pappa_egenes_qtls$variantID) |> 
  dplyr::select(Indiv, gt_GT_alleles, ID) |> 
  dplyr::rename(Donor = Indiv,
                variantID = ID)

# Normalized expression data for PAPPA
CTL_normQuant <- read_delim("data/rna/CTL_CPMadjTMM_invNorm.bed.gz") |> 
  filter(gene_id == pappa_egene_pbs$gene_id) |> 
  dplyr::select(-`#chr`, -start, -end, -length, -strand) |> 
  pivot_longer(cols = starts_with("AM"), names_to = "Donor", values_to = "expression") |> 
  mutate(Condition = "PBS")

FNF_normQuant <- read_delim("data/rna/FNF_CPMadjTMM_invNorm.bed.gz") |> 
  filter(gene_id == pappa_egene_fnf$gene_id)  |> 
  dplyr::select(-`#chr`, -start, -end, -length, -strand) |> 
  pivot_longer(cols = starts_with("AM"), names_to = "Donor", values_to = "expression") |> 
  mutate(Condition = "FN-f")

ALL_normQuant <- bind_rows(CTL_normQuant, FNF_normQuant)

# Join all data elements together
pappa_boxplot_data_qtl <- left_join(pappa_egenes_qtls, geno_data, by = "variantID") |> 
  left_join(ALL_normQuant, by = c("gene_id", "Donor"), relationship = "many-to-many")

# Get genotype orders
geno_order_qtls <- lapply(unique(pappa_boxplot_data_qtl$variantID),
                          determine_geno_order,
                          data = pappa_boxplot_data_qtl) |> 
  bind_rows()
# Join back geno orders
pappa_boxplot_data_qtl <- pappa_boxplot_data_qtl |> 
  left_join(geno_order_qtls, by = c("variantID", "gt_GT_alleles"))

# Put condition in proper order
pappa_boxplot_data_qtl$Condition <- factor(pappa_boxplot_data_qtl$Condition, 
                                           levels = c("PBS", "FN-f"))


# Create boxplots
pappa_boxplot_pbs_qtl <- create_eqtl_boxplot_horizontal(pappa_boxplot_data_qtl,
                                                        group = "PBS", 
                                                        stat_loc = "top",
                                                        rsID_loc = "top",
                                                        highlightAllele = NA,
                                                        condition_labs = "bottom")

pappa_boxplot_fnf_qtl <- create_eqtl_boxplot_horizontal(pappa_boxplot_data_qtl,
                                                        group = "FN-f", 
                                                        stat_loc = "top",
                                                        rsID_loc = "top",
                                                        highlightAllele = NA,
                                                        condition_labs = "bottom")
# Remove y axis title, labels, and ticks from FN-f plot to align with PBS
pappa_boxplot_fnf_qtl <- pappa_boxplot_fnf_qtl +
  theme(axis.line.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())


# Join PBS and FN-f boxplots
pappa_boxplots_qtl <- pappa_boxplot_pbs_qtl + pappa_boxplot_fnf_qtl +
  plot_annotation(theme = theme(panel.background = 
                                  element_rect(fill = "transparent", 
                                               color = "transparent"),
                                plot.background = 
                                  element_rect(fill = "transparent", 
                                               color = "transparent")))

save(pappa_boxplots_qtl, 
     file = "plots/pappaExample_Fig6_supp/pappa_boxplots_qtl.rda")


# PAPPA FN-f eQTL with GWAS ---------------------------------------------------------

pappa_colocs <- read_csv("data/colocalization/processed/PBS_FNF_allSignals_EUR_GWAS_250Kb_colocs.csv") |> 
  filter(eGene_name == "PAPPA")

pappa_pbs_rsid <- pappa_colocs |> filter(Condition == "PBS") |> pull(qtl_rsid)
pappa_fnf_rsid <- pappa_colocs |> filter(Condition == "FN-f") |> pull(qtl_rsid)


gwas_pappa <- read_csv(paste0("/work/users/n/e/nekramer/External/gwas/Boer_reprocessed/",
                              unique(pappa_colocs$OA), "/leads/EUR_", 
                              unique(pappa_colocs$OA), "_leads_ld_final.csv"),
                       col_types = "cdddccccddddddddddddddddddddddddddddddddddddddddddccdcccd") |> 
  filter(ldbuddy_rsID == pappa_pbs_rsid | ldbuddy_rsID == pappa_fnf_rsid) |> 
  dplyr::select(rsID, `CHR:hg38POS`, EA, NEA) |> 
  distinct() |> 
  mutate(variantID_v1 = paste0("chr", `CHR:hg38POS`, ":", EA, ":", NEA),
         variantID_v2 = paste0("chr", `CHR:hg38POS`, ":", NEA, ":", EA))

pappa_chrom <- unlist(str_split(gwas_pappa$`CHR:hg38POS`, ":"))[1]

# Check which variantID version is in our qtl dataset
chrom_fnf_nom <- fread(paste0("data/eqtl/qtl_nom/FNF_PEER_k22_genoPC_allSignals_nom1Mb_MAFs_chr", pappa_chrom, ".csv"), 
                       data.table = FALSE)
gwas_pappa_variantID <- c(gwas_pappa$variantID_v1, gwas_pappa$variantID_v2)[which(c(gwas_pappa$variantID_v1, gwas_pappa$variantID_v2) %in% chrom_fnf_nom$variantID)]

gwas_pappa <- gwas_pappa |> 
  pivot_longer(cols = contains("variantID"), 
               names_to = "variantID_version", 
               values_to = "variantID") |> 
  dplyr::select(-variantID_version) |> 
  filter(variantID == gwas_pappa_variantID)


pappa_egene <- read_csv("data/eqtl/FNF_PEER_k22_genoPC_cond1Mb_topSignals_rsID.csv") |>
  filter(gene_symbol == "PAPPA" & signal == 0)
pappa_egene_ld <- fread("data/eqtl/FNF_PEER_k22_genoPC_cond1Mb_topSignals_rsID_LD.csv",
                        data.table = FALSE) |>
  filter(gene_symbol == "PAPPA" & signal == 0) |>
  dplyr::select(ld_variantID, R2)

# GWAS
# Signal with LD buddies
gwas_pappa_data <- read_csv(paste0("/work/users/n/e/nekramer/External/gwas/Boer_reprocessed/",
                                   unique(pappa_colocs$OA), "/leads/EUR_",
                                   unique(pappa_colocs$OA), "_leads_ld_final.csv"),
                            col_types = "cdddccccddddddddddddddddddddddddddddddddddddddddddccdcccd") |> 
  filter(rsID == gwas_pappa$rsID) |> 
  filter(chrom == gsub("chr", "", pappa_egene[["gene_chr"]])) |> 
  mutate(ldbuddy_hg38pos = as.integer(str_extract(`ldbuddy_CHR:hg38POS`, "(?<=:)[0-9]+"))) |> 
  filter(ldbuddy_hg38pos >= signif((as.numeric(pappa_egene[["variant_start"]]) - 600000), digits = 5) 
         & ldbuddy_hg38pos <= signif((as.numeric(pappa_egene[["variant_start"]]) + 400000), digits = 5)) |> 
  dplyr::select(ldbuddy_rsID, `ldbuddy_CHR:hg38POS`, chrom, ldbuddy_hg38pos, p, ldbuddy_R2) |> 
  dplyr::rename(snp = ldbuddy_rsID, pos = ldbuddy_hg38pos)


# Get summary stats for other variants in region w/o R2 info
pappa_summary_stats <- fread(paste0("/work/users/n/e/nekramer/External/gwas/Boer_reprocessed/",
                                    unique(pappa_colocs$OA), "/summary_stats/", 
                                    unique(pappa_colocs$OA), "_", 
                                    pappa_egene[["gene_chr"]], ".csv"),
                             data.table = FALSE) |> 
  filter(hg38pos >= signif((as.numeric(pappa_egene[["variant_start"]]) - 600000), digits = 5) & 
           hg38pos <= signif((as.numeric(pappa_egene[["variant_start"]]) + 400000), digits = 5)) |>
  filter(!`CHR:hg38POS` %in% gwas_pappa_data$`ldbuddy_CHR:hg38POS`) |> 
  dplyr::select(`CHR:hg38POS`, chrom, hg38pos, p) |> 
  dplyr::rename(snp = `CHR:hg38POS`, pos = hg38pos) |> 
  mutate(ldbuddy_R2 = NA)

pappa_gwas_signal <- gwas_pappa_data |> 
  dplyr::select(-`ldbuddy_CHR:hg38POS`) |> 
  bind_rows(pappa_summary_stats) |> 
  mutate(chrom = paste0("chr", chrom)) |> 
  # Create group for LD group
  mutate(LDgrp = factor(cut(ldbuddy_R2, c(0, 0.2, 0.4, 0.6, 0.8, 1)),
                        levels = c("(0,0.2]", "(0.2,0.4]", "(0.4,0.6]",
                                   "(0.6,0.8]", "(0.8,1]", NA)))
pappa_gwas_signal$LDgrp <- addNA(pappa_gwas_signal$LDgrp)

# eQTL
pappa_fnf_signal <- fread(paste0("data/eqtl/qtl_nom/FNF_PEER_k22_genoPC_allSignals_nom1Mb_MAFs_",
                                 pappa_egene[["gene_chr"]], ".csv"), data.table = FALSE) |> 
  filter(gene_id == pappa_egene[["gene_id"]] & signal == 0) |> 
  mutate(rsID = case_when(variantID == pappa_egene[["variantID"]] ~ pappa_egene[["rsID"]],
                          variantID == gwas_pappa$variantID ~ gwas_pappa[["rsID"]])) |> 
  dplyr::rename(chrom = variant_chr,
                pos = variant_start,
                p = nom_pval,
                snp = rsID) |> 
  left_join(pappa_egene_ld, by = join_by(variantID == ld_variantID)) |> 
  mutate(LDgrp = cut(R2, c(0, 0.2, 0.4, 0.6, 0.8, 1)))


pappa_fnf_signal$LDgrp <- addNA(pappa_fnf_signal$LDgrp)


# plotgardener layout -----------------------------------------------------

pdf("plots/pappaExample_Fig6_supp/SupFig9.pdf", width = 4.55, height = 5.2)
pageCreate(width = 4.55, height = 5.2, showGuides = FALSE)
plotText("A", x = 0.1, y = 0.1, just = c("left", "top"), fontfamily = "Helvetica",
         fontsize = 14, fontface = "bold")
plotGG(pappa_boxplots_qtl, x = 0.25, y = 0.1, width = 4.3, height = 2)
plotText(label = "PBS lead variant", x = 1.5, y = 0.2, fontfamily = "Helvetica",
         fontsize = 8)
plotText(label = "FN-f lead variant", x = 3.5, y = 0.2, fontfamily = "Helvetica",
         fontsize = 8)


# FN-f eQTL signal
plotText("B", x = 0.1, y = 2, just = c("left", "top"), fontfamily = "Helvetica",
         fontsize = 14, fontface = "bold")

pappa_region <- pgParams(assembly = "hg38",
                         chrom = pappa_egene[["gene_chr"]],
                         chromstart = signif(as.numeric(pappa_egene[["variant_start"]]) - 600000, digits = 5),
                         chromend = signif(as.numeric(pappa_egene[["variant_start"]] + 400000), digits = 5),
                         x = 0.55,
                         width = 3.9,
                         fontfamily = "Helvetica")

qtl_pappa_ylim <- ceiling(max(-1*log10(pappa_signals$p))) + 1.5
gwas_pappa_ylim <- ceiling(max(-1*log10(pappa_gwas_signal$p))) + 1.5

# GWAS
gwas_pappa_plot <- plotManhattan(data = pappa_gwas_signal |> 
                                   arrange(desc(LDgrp)),
                                 params =  pappa_region,
                                 range = c(0, gwas_pappa_ylim),
                                 y = 3.15, height = 1, 
                                 sigline = FALSE,just = c("left", "bottom"),
                                 fill = colorby("LDgrp",
                                                palette = colorRampPalette(c("#262C74",
                                                                             "#98CDED",
                                                                             "#499A53",
                                                                             "#EEA741",
                                                                             "#DD3931",
                                                                             "grey"))),
                                 snpHighlights = data.frame(snp = c(pappa_fnf_rsid, gwas_pappa$rsID),
                                                            pch = c(24, 23),
                                                            cex = c(0.75, 0.75),
                                                            col = c("black", "black")))
annoYaxis(plot = gwas_pappa_plot, at = seq(0, gwas_pappa_ylim, 4), 
          axisLine = TRUE, fontsize = 6, gp = gpar(fontfamily = "Helvetica"))

plotText(
  label = "-log10(p-value)", x = 0.275, y = 2.6, rot = 90,
  fontsize = 6, just = "center",
  default.units = "inches", fontfamily = "Helvetica"
)
plotText(label = "Total Hip Replacement\nGWAS", x = 4.45, y = 2.1, 
         just = c("right", "top"), fontface = "bold",
         fontsize = 8, fontfamily = "Helvetica", lineheight = 0.8)

# LD legend
grid.points(x = 0.7, y = 2.25, default.units = "native", pch = 23, size = unit(0.5, "char"))
plotText(label = paste0(gwas_pappa$rsID, " (GWAS index)"),
         fontsize = 7, fontfamily = "Helvetica",
         just = "left", x = 0.8, y = 2.25)
plotLegend(legend = c("0.8 - 1.0",
                      "0.6 - 0.8",
                      "0.4 - 0.6",
                      "0.2 - 0.4",
                      "0.0 - 0.2"),
           fill = c("#DD3931", "#EEA741", "#499A53","#98CDED","#262C74"),
           x = 0.75, y = 2.45, width = 0.1, height = 0.4, border = FALSE, 
           fontsize = 6)
plotText(label = "r2", x = 1.375, y = 2.4, fontsize = 6)
# FN-f eQTL
fnf_pappa <- plotManhattan(data = pappa_fnf_signal,
                           params = pappa_region, 
                           range = c(0, qtl_pappa_ylim),
                           snpHighlights = data.frame(snp = c(gwas_pappa$rsID,
                                                              pappa_fnf_rsid),
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
                           y = 3.3, height = 1)
annoYaxis(plot = fnf_pappa, at = seq(0, qtl_pappa_ylim, 2), 
          axisLine = TRUE, fontsize = 6, gp = gpar(fontfamily = "Helvetica"))
plotText(
  label = "-log10(p-value)", x = 0.275, y = 3.85, rot = 90,
  fontsize = 6, just = "center",
  default.units = "inches", fontfamily = "Helvetica"
)
plotText(label = "eQTL", x = 4.45, y = 3.4, just = "right", fontface = "bold",
         fontsize = 8, fontfamily = "Helvetica")
plotText(label = paste0("PP4 = ", 
                        signif(pappa_colocs |> filter(Condition == "FN-f") |> pull(PP4_fnf), digits = 3)), 
         x = 4.45, y = 3.5, 
         just = "right", fontsize = 7, fontfamily = "Helvetica")


# LD legend
grid.points(x = 0.7, y = 3.5, default.units = "native", pch = 24, size = unit(0.45, "char"))
plotText(label = paste0(pappa_fnf_rsid, " (eQTL index)"),
         fontsize = 7, fontfamily = "Helvetica",
         just = "left", x = 0.8, y = 3.475)
plotLegend(legend = c("0.8 - 1.0",
                      "0.6 - 0.8",
                      "0.4 - 0.6",
                      "0.2 - 0.4",
                      "0.0 - 0.2"),
           fill = c("#DD3931", "#EEA741", "#499A53","#98CDED","#262C74"),
           x = 0.75, y = 3.7, width = 0.1, height = 0.4, border = FALSE, 
           fontsize = 6)
plotText(label = "r2", x = 1.375, y = 3.65, fontsize = 6)

pappa_gene_plot <- plotGenes(params = pappa_region, 
                             geneHighlights = data.frame("gene" = "PAPPA",
                                                         "color" = yl_gn_bu[6]),
                             y = 4.4, height = 0.5)
annoGenomeLabel(plot = pappa_gene_plot, y = 4.9,
                x = 0.55, fontsize = 8)

dev.off()
