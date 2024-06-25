library(tidyverse)
library(ggvenn)
library(vcfR)
library(ggtext)
library(patchwork)
library(plotgardener)
source("../plotting_utils.R")
source("../utils.R")


# Functions ---------------------------------------------------------------

create_eqtl_boxplot <- function(group, data){
  
  # Filter data for data from that grouping
  boxplot_data <- data |> 
    filter(grouping == group) |> 
    arrange(gt_order) |> 
    mutate(gt_GT_alleles = fct_inorder(gt_GT_alleles))
  
  # Grab rsID and eGene name for labeling
  rsid <- unique(boxplot_data$rsID)
  eGene_name <- unique(boxplot_data$gene_symbol)
  
  # Create dummy data to label betas and pvals 
  stat_data <- unique(boxplot_data[,c("beta_pbs", "beta_fnf", "nompval_pbs", "nompval_fnf")]) |> 
    pivot_longer(everything()) |> 
    separate_wider_delim(cols = "name", delim = "_", names = c("stat", "Condition")) |> 
    mutate(Condition = toupper(Condition)) |> 
    mutate(Condition = ifelse(Condition == "FNF", "FN-f", Condition)) |> 
    pivot_wider(names_from = stat, values_from = value) |> 
    mutate(Condition = factor(Condition, levels = c("PBS", "FN-f"))) |> 
    mutate(gt_GT_alleles = 3,
           beta_expression = 3.5,
           nom_expression = 3.35)
    
  
  eqtl_boxplot <- ggplot(boxplot_data, aes(x = gt_GT_alleles, y = expression)) +
    geom_hline(yintercept = 0, lty = 2, linewidth = 0.3) +
    geom_boxplot(aes(color = Condition, fill = Condition), outlier.shape = NA) +
    geom_jitter(position = position_jitter(width = 0.1), color = "grey25", size = 0.25) +
    # beta labels
    geom_text(data = stat_data, aes(y = beta_expression, label = paste0("Beta = ", signif(beta, digits = 3))),
              family = "Helvetica", vjust = 0, size = 3, hjust = 1) +
    geom_text(data = stat_data, aes(y = nom_expression, label = paste0("p-value = ", signif(nompval, digits = 3))),
              family = "Helvetica", vjust = 1, size = 3, hjust = 1) +
    scale_color_manual(values = c("#48617b", "#97723e")) +
    scale_fill_manual(values = c("#78A1Cd", "#FBBE67"))  + 
    theme_custom_general() + 
    xlab(rsid) +
    facet_wrap(vars(Condition), nrow = 2, strip.position = "left") +
    scale_y_continuous(name = paste0("**", eGene_name, "**", " normalized expression"),
                       limits = c(-3, 3.5), breaks = seq(-3, 3, 1)) +
    theme(legend.position = "none",
          axis.line.y = element_line(linewidth = 0.3),
          axis.ticks.y = element_line(linewidth = 0.3),
          strip.background = element_blank(),
          plot.background = element_rect(fill = "transparent", color = "transparent"),
          panel.background = element_rect(fill = "transparent", color = "transparent"),
          strip.placement = "outside",
          strip.text.y.left = element_blank(),
          strip.text = element_text(color = "black", size = 12),
          axis.text = element_text(color = "black"),
          text = element_text(family = "Helvetica", size = 10),
          axis.title.y = element_markdown(size = 9),
          plot.title = element_text(hjust = 0.5)) +
    coord_cartesian(clip = "off") +
    ggtitle(group)
  
  
  return(eqtl_boxplot)
}

# PBS, FN-f, and response Venn diagram ------------------------------------

## All eGenes
PBS_eGenes <- read_csv("data/eqtl/CTL_PEER_k20_genoPC_perm1Mb_sig_rsID.csv") 
FNF_eGenes <- read_csv("data/eqtl/FNF_PEER_k22_genoPC_perm1Mb_sig_rsID.csv") 

pbs_fnf_egenes <- tibble(values = unique(c(PBS_eGenes$gene_id, 
                                           FNF_eGenes$gene_id))) |> 
  mutate(PBS = values %in% PBS_eGenes$gene_id,
         FNF = values %in% FNF_eGenes$gene_id)


egene_venn <- ggplot(pbs_fnf_egenes, aes(A = PBS, B = FNF)) +
  geom_venn(set_names = c("PBS eGenes", "FN-f eGenes"), 
            fill_color = c(log2fcColors[["-"]], log2fcColors[["+"]]), 
            fill_alpha = 0.3,
            stroke_color = NA, auto_scale = TRUE, show_percentage = FALSE,
            text_size = 3, set_name_size = 3.25) +
  coord_fixed()  +
  theme(panel.background = element_rect(fill = "transparent", color = "transparent"),
        plot.background = element_rect(fill = "transparent", color = "transparent"),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank())
egene_venn <- venn_font(egene_venn, font = "Helvetica")

save(egene_venn, file = "plots/reQTLs_Fig3/egene_venn.rda")

## PBS:genotype interaction and specific pie chart
pbs_response_eGenes <- read_csv("data/reqtl/CTL_sig01_reQTLs_PEER_k20_genoPC.csv")
pbs_highconf_response_eGenes <- read_csv("data/reqtl/CTL_sig01_beta_donor_reQTLs_PEER_k20_genoPC.csv")

pbs_response_pie_data <- tibble(group = c("significant PBS:genotype interaction", 
                                          "high confidence PBS-specific"),
                                count = c(nrow(pbs_response_eGenes), nrow(pbs_highconf_response_eGenes))) |> 
  # Calculate percentages
  mutate(prop = count/nrow(pbs_response_eGenes)*100) |> 
  # Calculate position for text labels
  mutate(text_pos = cumsum(prop) - 0.5*prop) |> 
  mutate(count_label = ifelse(group == "high confidence PBS-specific", count, "")) |> 
  mutate(group_label = ifelse(group == "high confidence PBS-specific", 
                              "high\nconfidence\nPBS specific", ""))

pbs_response_all_col <- plotgardener:::makeTransparent(color = log2fcColors[["-"]], alpha = 0.6)
pbs_response_highconf_col <- log2fcColors[["-"]]

pbs_response_pie_chart <- ggplot(pbs_response_pie_data, aes(x = "", y = prop, fill = group)) +
  geom_bar(stat = "identity", width = 1) +
  scale_fill_manual(values = c(pbs_response_highconf_col, pbs_response_all_col)) +
  coord_polar("y", start = 2*pi/3, clip = "off") +
  geom_text(aes(y = text_pos, label = count_label), family = "Helvetica", 
            size = 3, fontface = "bold") +
  geom_text(aes(x = 1.55, label = group_label), family = "Helvetica",
            size = 3, position = position_stack(vjust = 0.5), hjust = 0) +
  theme(legend.position = "none",
        panel.background = element_rect(fill = 'transparent', color = 'transparent'),
        plot.background = element_rect(fill = 'transparent', color = 'transparent'),
        plot.title = element_markdown(size = 9, family = "Helvetica", hjust = 0.5, 
                                      margin = margin(b = -10)),
        axis.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()) +
  ggtitle(label = paste0("**", nrow(pbs_response_eGenes), "**"))

save(pbs_response_pie_chart, file = "plots/reQTLs_Fig3/pbs_response_pie_chart.rda")

## FNF:genotype interaction and response pie chart
fnf_response_eGenes <- read_csv("data/reqtl/FNF_sig01_reQTLs_PEER_k22_genoPC.csv")
fnf_highconf_response_eGenes <- read_csv("data/reqtl/FNF_sig01_beta_donor_reQTLs_PEER_k22_genoPC.csv")

fnf_response_pie_data <- tibble(group = c("FN-f response", 
                                          "high confidence FN-f response"),
                                count = c(nrow(fnf_response_eGenes), 
                                          nrow(fnf_highconf_response_eGenes))) |> 
  # Calculate percentages
  mutate(prop = count/nrow(fnf_response_eGenes)*100) |> 
  # Calculate position for text labels
  mutate(text_pos = cumsum(prop) - 0.5*prop) |> 
  mutate(count_label = ifelse(group == "high confidence FN-f response", count, "")) |> 
  mutate(group_label = ifelse(group == "high confidence FN-f response", 
                              "high\nconfidence\nFN-f response", ""))
fnf_response_pie_data$group <- factor(fnf_response_pie_data$group, 
                                      levels = c("high confidence FN-f response",
                                                 "FN-f response"))

fnf_response_all_col <- plotgardener:::makeTransparent(color = log2fcColors[["+"]], alpha = 0.6)
fnf_response_highconf_col <- log2fcColors[["+"]]

fnf_response_pie_chart <- ggplot(fnf_response_pie_data, aes(x = "", y = prop, fill = group)) +
  geom_bar(stat = "identity", width = 1) +
  scale_fill_manual(values = c(fnf_response_highconf_col, fnf_response_all_col)) +
  coord_polar("y", start = 7*pi/4, clip = "off") + 
  geom_text(aes(y = text_pos, label = count_label), family = "Helvetica",
            size = 3, fontface = "bold") +
  geom_text(aes(x = 1.55, label = group_label), family = "Helvetica",
            size = 3, position = position_stack(vjust = 0.5), hjust = 1) +
  theme(legend.position = "none",
        panel.background = element_rect(fill = 'transparent', color = 'transparent'),
        plot.background = element_rect(fill = 'transparent', color = 'transparent'),
        plot.title = element_markdown(size = 9, family = "Helvetica", hjust = 0.5, 
                                      margin = margin(b = -10)),
        axis.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())  +
  ggtitle(label = paste0("**", nrow(fnf_response_eGenes), "**"))

save(fnf_response_pie_chart, file = "plots/reQTLs_Fig3/fnf_response_pie_chart.rda")

# boxplots - shared (GPX7), PBS-specific (MPZL2), and FN-f-specific (DIO2)-----------------

pbs_egenes <- read_csv("data/eqtl/CTL_PEER_k20_genoPC_perm1Mb_sig_rsID.csv",
                       col_select = c("rsID", "gene_id", "variantID", "gene_symbol")) |>
  filter(gene_symbol %in% c("GPX7", "MPZL2")) |> 
  mutate(grouping = c("shared", "PBS-specific"),
         signal = 0)
fnf_egenes <- read_csv("data/eqtl/FNF_PEER_k22_genoPC_perm1Mb_sig_rsID.csv",
                       col_select = c("rsID", "gene_id", "variantID", "gene_symbol")) |> 
  filter(gene_symbol == "DIO2") |> 
  mutate(grouping = c("FN-f-specific"),
         signal = 0)

boxplot_egenes <- bind_rows(pbs_egenes, fnf_egenes)
# Get minor allele
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
  left_join(ALL_normQuant, by = c("gene_id", "Donor"))
  
# Get orders of genotypes  
geno_orders <- lapply(unique(egene_boxplot_data$variantID), 
                        determine_geno_order, 
                        data = egene_boxplot_data) |> 
  bind_rows()
  
egene_boxplot_data <- egene_boxplot_data |> 
  left_join(geno_orders, by = c("variantID", "gt_GT_alleles"))
egene_boxplot_data$Condition <- factor(egene_boxplot_data$Condition, levels = c("PBS", "FN-f"))
egene_boxplot_data$grouping <- factor(egene_boxplot_data$grouping, levels = c("shared", "PBS-specific", "FN-f-specific"))

# Create boxplots
shared_boxplot <- create_eqtl_boxplot(group = "shared", egene_boxplot_data)
pbs_boxplot <- create_eqtl_boxplot(group = "PBS-specific", egene_boxplot_data)
fnf_boxplot <- create_eqtl_boxplot(group = "FN-f-specific", egene_boxplot_data)

reqtl_boxplots <- shared_boxplot + pbs_boxplot + fnf_boxplot +
  plot_annotation(theme = theme(panel.background = 
                                  element_rect(fill = "transparent", 
                                               color = "transparent"),
                                plot.background = 
                                  element_rect(fill = "transparent", 
                                               color = "transparent")))

save(reqtl_boxplots, file = "plots/reQTLs_Fig3/reqtl_boxplots.rda")

# Data for locus zooms ----------------------------------------------------

pbs_egene <- read_csv("data/eqtl/CTL_PEER_k20_genoPC_perm1Mb_sig_rsID.csv") |>
  filter(gene_symbol == "MRPL50")

pbs_egene_ld <- fread("data/eqtl/CTL_PEER_k20_genoPC_perm1Mb_sig_rsID_LD_rsID.csv",
                      data.table = FALSE) |> 
  filter(gene_symbol == "MRPL50") |> 
  dplyr::select(ld_variantID, R2)

mrpl50_signals <- bind_rows(fread(paste0("data/eqtl/CTL_PEER_k20_genoPC_nom1Mb_MAFs_", 
                                         pbs_egene[["gene_chr"]], ".csv"),
                                  data.table = FALSE) |> 
                              filter(gene_id == pbs_egene[["gene_id"]]) |> 
                              mutate(Condition = "PBS"),
                            fread(paste0("data/eqtl/FNF_PEER_k22_genoPC_nom1Mb_MAFs_", 
                                         pbs_egene[["gene_chr"]], ".csv"),
                                  data.table = FALSE) |> 
                              filter(gene_id == pbs_egene[["gene_id"]]) |> 
                              mutate(Condition = "FN-f")) |>
  mutate(rsID = ifelse(variantID == pbs_egene[["variantID"]], pbs_egene[["rsID"]], NA)) |> 
  dplyr::rename(chrom = variant_chr,
                pos = variant_start,
                p = nom_pval,
                snp = rsID) |> 
  left_join(pbs_egene_ld, by = join_by(variantID == ld_variantID)) |> 
  mutate(LDgrp = cut(R2, c(0, 0.2, 0.4, 0.6, 0.8, 1)))

mrpl50_signals$LDgrp <- addNA(mrpl50_signals$LDgrp)

fnf_egene <-  read_csv("data/eqtl/FNF_PEER_k22_genoPC_perm1Mb_sig_rsID.csv") |>
  filter(gene_symbol == "SAP30BP")

fnf_egene_ld <- fread("data/eqtl/FNF_PEER_k22_genoPC_perm1Mb_sig_rsID_LD_rsID.csv",
                      data.table = FALSE) |> 
  filter(gene_symbol == "SAP30BP") |> 
  dplyr::select(ld_variantID, R2)

sap30bp_signals <- bind_rows(fread(paste0("data/eqtl/CTL_PEER_k20_genoPC_nom1Mb_MAFs_", 
                                          fnf_egene[["gene_chr"]], ".csv"),
                                   data.table = FALSE) |> 
                               filter(gene_id == fnf_egene[["gene_id"]]) |> 
                               mutate(Condition = "PBS"),
                             fread(paste0("data/eqtl/FNF_PEER_k22_genoPC_nom1Mb_MAFs_", 
                                          fnf_egene[["gene_chr"]], ".csv"),
                                   data.table = FALSE) |> 
                               filter(gene_id == fnf_egene[["gene_id"]]) |> 
                               mutate(Condition = "FN-f")) |>
  mutate(rsID = ifelse(variantID == fnf_egene[["variantID"]], fnf_egene[["rsID"]], NA)) |> 
  dplyr::rename(chrom = variant_chr,
                pos = variant_start,
                p = nom_pval,
                snp = rsID) |> 
  left_join(fnf_egene_ld, by = join_by(variantID == ld_variantID)) |> 
  mutate(LDgrp = cut(R2, c(0, 0.2, 0.4, 0.6, 0.8, 1))) 

sap30bp_signals$LDgrp <- addNA(sap30bp_signals$LDgrp)


# plotgardener layout -----------------------------------------------------

pdf("plots/reQTLs_Fig3/Fig3.pdf", width = 9.25, height = 7.75)

pageCreate(width = 9.25, height = 7.75, showGuides = FALSE)

### A - Venn diagrams/pie charts of response eQTLs
plotText("A", x = 0.1, y = 0.1, just = c("left", "top"), fontfamily = "Helvetica",
         fontsize = 14, fontface = "bold")

plotGG(egene_venn, x = -0.4, y = -0.4, width = 4, height = 3.25)

plotSegments(x0 = 0.35, x1 = 0.35, y0 = 1.25, y1 = 2.05, lty = 2)
plotSegments(x0 = 0.35, x1 = 0.5, y0 = 1.25, y1 = 1.25, lty = 2)

plotSegments(x0 = 2.9, x1 = 2.9, y0 = 1.25, y1 = 2.05, lty = 2)
plotSegments(x0 = 2.9, x1 = 2.75, y0 = 1.25, y1 = 1.25, lty = 2)

plotText(label = "significant genetic and condition interaction effect",
         fontfamily = "Helvetica", fontsize = 9, x = 1.6, y = 2.25)

plotGG(pbs_response_pie_chart, x = -0.1, y = 2.523, width = 1.25*(696/856), 
       height = 1.25*(696/856))
plotGG(fnf_response_pie_chart, x = 2.1, y = 2.4, width = 1.25, height = 1.25)

## B
plotText("B", x = 3.2, y = 0.1, just = c("left", "top"), fontfamily = "Helvetica",
         fontsize = 14, fontface = "bold")
plotGG(plot = reqtl_boxplots, x = 3.3, y = -0.1, width = 6.1, height = 3.75)
plotText(label = "PBS", x = 3.6, y = 0.75, just = "right", fontsize = 9,
         fontfamily = "Helvetica")
plotText(label = "FN-f", x = 3.6, y = 2.75, just = "right", fontsize = 9,
         fontfamily = "Helvetica")

## C PBS-specific signal
# MRPL50
plotText("C", x = 0.1, y = 3.65, just = "left", fontfamily = "Helvetica",
         fontsize = 14, fontface = "bold")
mrpl50_region <- pgParams(assembly = "hg38", gene = "MRPL50", geneBuffer = 500000, 
                          x = 0.55, width = 3.75)
mrpl50_ylim <- ceiling(max(log10(mrpl50_signals$p)*-1)) + 2

pbs_mrpl50 <- plotManhattan(data = mrpl50_signals |> filter(Condition == "PBS"),
                            params = mrpl50_region, 
                            range = c(0, mrpl50_ylim),
                            snpHighlights = data.frame(snp = pbs_egene[["rsID"]],
                                                       pch = 23,
                                                       cex = 0.5),
                            fill = colorby("LDgrp",
                                           palette = colorRampPalette(c("#262C74",
                                                                        "#98CDED",
                                                                        "#499A53",
                                                                        "#EEA741",
                                                                        "#DD3931",
                                                                        "grey"))),
                            y = 3.65, height = 1.5)

annoYaxis(plot = pbs_mrpl50, at = seq(0, mrpl50_ylim, 2), 
          axisLine = TRUE, fontsize = 8)
plotText(
  label = "-log10(p-value)", x = 0.24, y = 4.4, rot = 90,
  fontsize = 8, just = "center",
  default.units = "inches", fontfamily = "Helvetica"
)
plotText(label = "PBS", x = 0.65, y = 3.65, just = c("left", "top"),
         fontfamily = "Helvetica", fontsize = 11)

plotText(label = pbs_egene[["rsID"]], x = 1.675, y = 3.9, just = "left",
         fontfamily = "Helvetica", fontsize = 9)

fnf_mrpl50 <- plotManhattan(data = mrpl50_signals |> filter(Condition == "FN-f"),
                            params = mrpl50_region, 
                            range = c(0, mrpl50_ylim),
                            snpHighlights = data.frame(snp =  pbs_egene[["rsID"]],
                                                       pch = 23,
                                                       cex = 0.5),
                            fill = colorby("LDgrp",
                                           palette = colorRampPalette(c("#262C74",
                                                                        "#98CDED",
                                                                        "#499A53",
                                                                        "#EEA741",
                                                                        "#DD3931",
                                                                        "grey"))),
                            y = 5.3, height = 1.5)

annoYaxis(plot = fnf_mrpl50, at = seq(0, mrpl50_ylim, 2), 
          axisLine = TRUE, fontsize = 8)

plotText(
  label = "-log10(p-value)", x = 0.24, y = 6.1, rot = 90,
  fontsize = 8, just = "center",
  default.units = "inches", fontfamily = "Helvetica"
)
plotText(label = "FN-f", x = 0.65, y = 5.3, just = c("left", "top"),
         fontfamily = "Helvetica", fontsize = 11)

mrpl50_genes <- plotGenes(params = mrpl50_region, y = 6.9,
                   height = 0.6, geneHighlights = data.frame("gene" = "MRPL50",
                                                             "color" = "#37a7db"))
annoGenomeLabel(plot = mrpl50_genes, params = mrpl50_region,
                y = 7.5)

## D FN-f specific signal (backup CDC42EP3)
# SAP30BP
plotText("D", x = 04.5, y = 3.65, just = "left", fontfamily = "Helvetica",
         fontsize = 14, fontface = "bold")
sap30bp_region <- pgParams(assembly = "hg38", gene = "SAP30BP", geneBuffer = 500000, 
                          x = 5, width = 3.75)
sap30bp_ylim <- ceiling(max(log10(sap30bp_signals$p)*-1)) + 2

pbs_sap30bp <- plotManhattan(data = sap30bp_signals |> filter(Condition == "PBS"),
                            params = sap30bp_region, 
                            range = c(0, sap30bp_ylim),
                            snpHighlights = data.frame(snp = fnf_egene[["rsID"]],
                                                       pch = 23,
                                                       cex = 0.5),
                            fill = colorby("LDgrp",
                                           palette = colorRampPalette(c("#262C74",
                                                                        "#98CDED",
                                                                        "#499A53",
                                                                        "#EEA741",
                                                                        "#DD3931",
                                                                        "grey"))),
                            y = 3.65, height = 1.5)

annoYaxis(plot = pbs_sap30bp, at = seq(0, sap30bp_ylim, 2), 
          axisLine = TRUE, fontsize = 8)
plotText(
  label = "-log10(p-value)", x = 4.675, y = 4.4, rot = 90,
  fontsize = 8, just = "center",
  default.units = "inches", fontfamily = "Helvetica"
)
plotText(label = "PBS", x = 5.1, y = 3.65, just = c("left", "top"),
         fontfamily = "Helvetica", fontsize = 11)



fnf_sap30bp <- plotManhattan(data = sap30bp_signals |> filter(Condition == "FN-f"),
                            params = sap30bp_region, 
                            range = c(0, sap30bp_ylim),
                            snpHighlights = data.frame(snp = fnf_egene[["rsID"]],
                                                       pch = 23,
                                                       cex = 0.5),
                            fill = colorby("LDgrp",
                                           palette = colorRampPalette(c("#262C74",
                                                                        "#98CDED",
                                                                        "#499A53",
                                                                        "#EEA741",
                                                                        "#DD3931",
                                                                        "grey"))),
                            y = 5.3, height = 1.5)

annoYaxis(plot = fnf_sap30bp, at = seq(0, sap30bp_ylim, 2), 
          axisLine = TRUE, fontsize = 8)

plotText(
  label = "-log10(p-value)", x = 4.675, y = 6.1, rot = 90,
  fontsize = 8, just = "center",
  default.units = "inches", fontfamily = "Helvetica"
)
plotText(label = "FN-f", x = 5.1, y = 5.3, just = c("left", "top"),
         fontfamily = "Helvetica", fontsize = 11)

plotText(label = fnf_egene[["rsID"]], x = 6.7, y = 5.5, just = "left",
         fontfamily = "Helvetica", fontsize = 9)

sap30bp_genes <- plotGenes(params = sap30bp_region, y = 6.9,
                          height = 0.6, geneHighlights = data.frame("gene" = "SAP30BP",
                                                                    "color" = "#37a7db"))
annoGenomeLabel(plot = sap30bp_genes, params = sap30bp_region,
                y = 7.5)

plotLegend(legend = c("0.8 - 1.0",
                      "0.6 - 0.8",
                      "0.4 - 0.6",
                      "0.2 - 0.4",
                      "0.0 - 0.2"),
           fill = c("#DD3931", "#EEA741", "#499A53","#98CDED","#262C74"),
           x = 7.75, y = 3.75, width = 0.1, height = 0.5, border = FALSE, 
           fontsize = 8)
plotText("r2", x = 8.3, y = 3.675, fontsize = 8, fontfamily = "Helvetica")

dev.off()
