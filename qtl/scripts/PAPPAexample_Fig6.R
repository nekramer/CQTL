library(tidyverse)
library(patchwork)
library(plotgardener)
library(mariner)
library(ggtext)
source("../utils.R")
source("../plotting_utils.R")

# Coloc info 

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
chrom_pbs_nom <- fread(paste0("data/eqtl/qtl_nom/CTL_PEER_k20_genoPC_allSignals_nom1Mb_MAFs_chr", pappa_chrom, ".csv"), 
                       data.table = FALSE)
gwas_pappa_variantID <- c(gwas_pappa$variantID_v1, gwas_pappa$variantID_v2)[which(c(gwas_pappa$variantID_v1, gwas_pappa$variantID_v2) %in% chrom_pbs_nom$variantID)]

gwas_pappa <- gwas_pappa |> 
  pivot_longer(cols = contains("variantID"), 
               names_to = "variantID_version", 
               values_to = "variantID") |> 
  dplyr::select(-variantID_version) |> 
  filter(variantID == gwas_pappa_variantID)



# eQTL boxplots of PAPPA GWAS variant -------------------------------------


pappa_egene_gwas <- read_csv("data/eqtl/CTL_PEER_k20_genoPC_cond1Mb_topSignals_rsID_signalRanges.csv",
                       col_select = c("gene_id", "gene_symbol", "signal", 
                                      "gene_chr", "min_signal_pos", "max_signal_pos")) |>
  filter(gene_symbol == "PAPPA") |> 
  mutate(variantID = gwas_pappa$variantID) |> 
  mutate(rsID = gwas_pappa$rsID)

# Get minor allele for correct genotype ordering
pappa_egene_gwas$MA <- get_minor_allele(varID = gwas_pappa$variantID, 
                                   data = pappa_egene_gwas)

# Get PBS and FNF betas and nominal p-vals
pappa_egene_gwas <- pappa_egene_gwas |> 
  bind_cols(get_pvals_betas(varID = gwas_pappa$variantID, data = pappa_egene_gwas))

# Read VCF from pappa_chrom
# Geno data from VCF file containing GWAS, PBS, and FN-f variants
vcf <- vcfR2tidy(read.vcfR(paste0("data/vcf/CQTL_COA_01_GDA8_COA2_01_COA3_01_GDA8_COA4_COA5_COA6_COA7_qtl_chr",
                                  pappa_chrom,".vcf.gz"), 
                           verbose = FALSE))
variants <- vcf$fix |> 
  dplyr::select(POS, ID, REF, ALT)
geno_data <- vcf$gt |> 
  left_join(variants, by = "POS", relationship = "many-to-many") |> 
  filter(ID == gwas_pappa_variantID) |> 
  dplyr::select(Indiv, gt_GT_alleles, ID) |> 
  dplyr::rename(Donor = Indiv,
                variantID = ID)

# Normalized expression data for PAPPA
CTL_normQuant <- read_delim("data/rna/CTL_CPMadjTMM_invNorm.bed.gz") |> 
  filter(gene_id == pappa_egene_gwas$gene_id) |> 
  dplyr::select(-`#chr`, -start, -end, -length, -strand) |> 
  pivot_longer(cols = starts_with("AM"), names_to = "Donor", values_to = "expression") |> 
  mutate(Condition = "PBS")

FNF_normQuant <- read_delim("data/rna/FNF_CPMadjTMM_invNorm.bed.gz") |> 
  filter(gene_id == pappa_egene_gwas$gene_id)  |> 
  dplyr::select(-`#chr`, -start, -end, -length, -strand) |> 
  pivot_longer(cols = starts_with("AM"), names_to = "Donor", values_to = "expression") |> 
  mutate(Condition = "FN-f")

ALL_normQuant <- bind_rows(CTL_normQuant, FNF_normQuant)

# Join all data elements together
pappa_boxplot_data_gwas <- left_join(pappa_egene_gwas, geno_data, by = "variantID") |> 
  left_join(ALL_normQuant, by = c("gene_id", "Donor"), relationship = "many-to-many")

# Get genotype order 
geno_order_gwas <- determine_geno_order(varID = gwas_pappa$variantID, 
                     data = pappa_boxplot_data_gwas)

pappa_boxplot_data_gwas <- pappa_boxplot_data_gwas |> 
  left_join(geno_order_gwas, by = c("variantID", "gt_GT_alleles"))
pappa_boxplot_data_gwas$Condition <- factor(pappa_boxplot_data_gwas$Condition, 
                                            levels = c("PBS", "FN-f"))

# Create boxplots
pappa_boxplot_gwas <- create_eqtl_boxplot_horizontal(pappa_boxplot_data_gwas,
                                                stat_loc = "top",
                                                highlightAllele = "minor",
                                                condition_labs = "bottom")
# Adjust parameters to fit this specific figure
pappa_boxplot_gwas$theme$axis.title.y$size <- 6
pappa_boxplot_gwas$theme$axis.text.x$margin <- margin(b = -5)
pappa_boxplot_gwas$layers[[6]]$data$expression <- c(-3.5, -3.5)

pappa_boxplot_gwas <- pappa_boxplot_gwas + 
  scale_y_continuous(limits = c(-3.5, 3.5), breaks = seq(-3, 3, 1), 
                     name = paste0("**", unique(pappa_boxplot_gwas$data$gene_symbol), "**", " normalized expression"))

save(pappa_boxplot_gwas, file = "plots/pappaExample_Fig6/pappa_boxplot_gwas.rda")


# PAPPA FN-f boxplots -----------------------------------------------------

# Get counts for FN-f
load("../DE/data/condition_de/differential_expression_dds.rda")
pappa_fnf <- get_gene_condition_Counts(gene = pappa_egene_gwas[["gene_id"]],
                                       dds = dds) |> 
  mutate(Condition = ifelse(Condition == "CTL", "PBS", "FN-f")) |> 
  mutate(Condition = factor(Condition, levels = c("PBS", "FN-f")))

## DE stats
pappa_de_l2fc <- read_csv("../DE/data/condition_de/de_genes_results.csv") |> 
  filter(symbol == "PAPPA") |> 
  pull(log2FoldChange)


pappa_fnf_boxplot <- ggplot(pappa_fnf, aes(x = Condition, 
                      y = log2(count),
                      fill = Condition)) +
  geom_boxplot(outlier.shape = NA,
               linewidth = 0.25) +
  annotate(geom = "segment", x = 1, xend = 2,
           y = 16.95, yend = 16.95, linewidth = 0.25) +
  annotate(geom = "richtext", label = paste0("log~2~FC = ", signif(pappa_de_l2fc, digits = 3)),
           x = 1.5, y = Inf, family = "Helvetica", size = 2, fill = NA, label.color = NA) +
  scale_fill_manual(values = c(log2fcColors[["-"]], log2fcColors[["+"]])) +
  scale_y_continuous( name = "**PAPPA** log~2~(normalized counts)",
                     limits = c(10, 17), breaks = seq(10, 17, 2)) +
  coord_cartesian(clip = "off") +
  theme(strip.placement = "outside",
        axis.line = element_line(linewidth = 0.25),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(color = "black", linewidth = 0.25),
        axis.ticks.length.y = unit(-0.1, "cm"),
        axis.title.x = element_blank(),
        axis.title.y = element_markdown(size = 6, family= "Helvetica",
                                        margin = margin(r = -0.5)),
        text = element_text(family = "Helvetica"),
        axis.text.y = element_text(color = "black", size = 6),
        axis.text.x = element_text(color = "black", size = 6, margin = margin(b = 0)),
        strip.background = element_blank(),
        strip.text.x.top = element_text(size = 8, margin = margin(b = 5)),
        panel.background = element_rect(fill = "transparent", color = "transparent"),
        plot.background = element_rect(fill = "transparent", color = "transparent"),
        panel.grid = element_blank(),
        legend.position = "none",
        panel.spacing.y = unit(0.5, "cm"))

save(pappa_fnf_boxplot, file = "plots/pappaExample_Fig6/pappa_fnf_boxplot.rda")

# PAPPA age boxplots ------------------------------------------------------

# Get counts for age de genes in PBS and FN-f
load("../DE/data/age_de/dds_age_ctl_lrt.rda")
load("../DE/data/age_de/dds_age_fnf_lrt.rda")

pappa_age_pbs <- get_gene_age_Counts(gene = pappa_egene_gwas[["gene_id"]], 
                                       dds = dds_age_ctl_lrt) |> 
  mutate(condition = "PBS")

pappa_age_fnf <- get_gene_age_Counts(gene = pappa_egene_gwas[["gene_id"]], 
                                       dds = dds_age_fnf_lrt) |> 
  mutate(condition = "FN-f")

pappa_age <- bind_rows(pappa_age_pbs, pappa_age_fnf) |> 
  mutate(condition = factor(condition, levels = c("PBS", "FN-f"))) |> 
  # Split into age groups
  mutate(Age_group = ifelse(Age < 50, "<50", ifelse(Age < 65, "50-65", "65+"))) |> 
  mutate(Age_group = factor(Age_group, levels = c("<50", "50-65", "65+")))


ctl_cluster_pval05 <- read_csv("../DE/data/age_de/ctl_age_pval05clusters.csv") |> 
  distinct(gene_id, .keep_all = TRUE) |> 
  dplyr::select(gene_id, symbol, cluster) |> 
  dplyr::rename(ctl_cluster = cluster)

fnf_cluster_pval05 <- read_csv("../DE/data/age_de/fnf_age_pval05clusters.csv") |> 
  distinct(gene_id, .keep_all = TRUE) |> 
  dplyr::select(gene_id, symbol, cluster) |> 
  dplyr::rename(fnf_cluster = cluster)

pappa_sig_age_clusters <- full_join(ctl_cluster_pval05, fnf_cluster_pval05, 
                                    by = c("gene_id", "symbol")) |> 
  filter(symbol == "PAPPA") |> 
  pivot_longer(cols = ends_with("cluster"), 
               names_to = "condition", 
               values_to = "cluster") |> 
  filter(!is.na(cluster)) |> 
  mutate(condition = gsub("_cluster", "", condition)) |> 
  mutate(condition = ifelse(condition == "ctl", "PBS", "FN-f")) |> 
  mutate(condition = factor(condition, levels = c("PBS", "FN-f"))) |> 
  mutate(cluster = as.character(cluster)) |> 
  mutate(maxCount = max(log2(pappa_age$count))) |> 
  group_by(gene_id, symbol, condition, cluster, maxCount) |> 
  tidyr::expand(Age_group = c("<50", "50-65", "65+")) |> 
  ungroup()

pappa_age_boxplot <- ggplot(pappa_age, aes(x = Age_group, y = log2(count))) +
  geom_boxplot(outlier.shape = NA,
               linewidth = 0.25, fill = ageClusterColors[["+"]]) +
  geom_line(data = pappa_sig_age_clusters,
            aes(x = Age_group, y = maxCount + 0.5, group = condition), inherit.aes = FALSE,
            linewidth = 0.25) +
  geom_star(data = pappa_sig_age_clusters,
            aes(x = 2, y = Inf), fill = ageClusterColors[["+"]],
            size = 1.5, starstroke = 0.25) +
  facet_wrap(vars(condition), scales = "free") +
  scale_x_discrete(name = "Age group") +
  scale_y_continuous(breaks = scales::breaks_pretty(3), 
                     name = "**PAPPA** log~2~(normalized counts)") +
  coord_cartesian(clip = "off") +
  theme(strip.placement = "outside",
        axis.line = element_line(linewidth = 0.25),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(color = "black", linewidth = 0.25),
        axis.ticks.length.y = unit(-0.1, "cm"),
        axis.title.x = element_text(size = 6, family = "Helvetica", margin = margin(t = -2)),
        axis.title.y = element_markdown(size = 6, family= "Helvetica",
                                        margin = margin(r = -0.5)),
        text = element_text(family = "Helvetica"),
        axis.text.y = element_text(color = "black", size = 6),
        axis.text.x = element_text(color = "black", size = 6, margin = margin(b = 0)),
        strip.background = element_blank(),
        strip.text.x.top = element_text(size = 8, margin = margin(b = 5)),
        panel.background = element_rect(fill = "transparent", color = "transparent"),
        plot.background = element_rect(fill = "transparent", color = "transparent"),
        panel.grid = element_blank(),
        legend.position = "none",
        panel.spacing.y = unit(0.5, "cm"))

save(pappa_age_boxplot, file = "plots/pappaExample_Fig6/pappa_age_boxplot.rda")


# Data for PAPPA GWAS and QTL locus zooms ---------------------------------

pappa_egene <- read_csv("data/eqtl/CTL_PEER_k20_genoPC_cond1Mb_topSignals_rsID.csv") |>
  filter(gene_symbol == "PAPPA" & signal == 0)
pappa_egene_ld <- fread("data/eqtl/CTL_PEER_k20_genoPC_cond1Mb_topSignals_rsID_LD.csv",
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
pappa_signals <- bind_rows(fread(paste0("data/eqtl/qtl_nom/CTL_PEER_k20_genoPC_allSignals_nom1Mb_MAFs_", 
                                        pappa_egene[["gene_chr"]], ".csv"),
                                 data.table = FALSE) |> 
                             filter(gene_id == pappa_egene[["gene_id"]] & signal == 0) |> 
                             mutate(Condition = "PBS"),
                           fread(paste0("data/eqtl/qtl_nom/FNF_PEER_k22_genoPC_allSignals_nom1Mb_MAFs_", 
                                        pappa_egene[["gene_chr"]], ".csv"),
                                 data.table = FALSE) |> 
                             filter(gene_id == pappa_egene[["gene_id"]] & signal == 0) |> 
                             mutate(Condition = "FN-f")) |>
  mutate(rsID = ifelse(variantID == pappa_egene[["variantID"]], pappa_egene[["rsID"]],
                       ifelse(variantID == pappa_egene[["variantID"]], pappa_egene[["rsID"]],
                              ifelse(variantID == gwas_pappa$variantID, gwas_pappa$rsID, NA)))) |> 
  dplyr::rename(chrom = variant_chr,
                pos = variant_start,
                p = nom_pval,
                snp = rsID) |> 
  left_join(pappa_egene_ld, by = join_by(variantID == ld_variantID)) |> 
  mutate(LDgrp = cut(R2, c(0, 0.2, 0.4, 0.6, 0.8, 1)))
pappa_signals$LDgrp <- addNA(pappa_signals$LDgrp)

# PAPPA loops
load("data/hic/CQTL_10kb_static_loops.rda")
CQTL_30kb_static_loops <- snapToBins(CQTL_10kb_static_loops, 30000)


ensembl_txdb <- makeTxDbFromEnsembl(organism = "Homo sapiens", release = 111)
seqlevelsStyle(ensembl_txdb) <- "UCSC"

pappa_transcripts <- AnnotationDbi::select(ensembl_txdb, 
                                             keys = pappa_egene_gwas[["gene_id"]],
                                             keytype = "GENEID",
                                             columns = "TXID")

pappa_promoters <- promoters(ensembl_txdb) |> 
  filter(tx_id %in% pappa_transcripts$TXID)


pappa_loop_overlaps <- linkOverlaps(CQTL_30kb_static_loops, 
                                      subject1 = GRanges(seqnames = pappa_egene_gwas[["gene_chr"]],
                                                         ranges = IRanges(start = as.numeric(pappa_egene_gwas[["min_signal_pos"]]),
                                                                          end = as.numeric(pappa_egene_gwas[["max_signal_pos"]]))), 
                                      subject2 = pappa_promoters)

pappa_loop_10kb <- CQTL_10kb_static_loops[unique(pappa_loop_overlaps$query)]
pappa_loop_30kb <- CQTL_30kb_static_loops[unique(pappa_loop_overlaps$query)]

# Assemble figure with plotgardener ----------------------------------------------

pdf(file = "plots/pappaExample_Fig6/Fig6.pdf", width = 4.5, height = 8.5)
pageCreate(width = 4.5, height = 8.5, showGuides = FALSE)

## A - PAPPA FN-f boxplot
plotText("A", x = 0.1, y = 0.1, just = c("left", "top"), fontfamily = "Helvetica",
         fontsize = 14, fontface = "bold")
plotGG(pappa_fnf_boxplot, x = 0.3, y = 0.1, width = 1.5, height = 1.5)


## B - PAPPA age expression boxplot
plotText("B", x = 1.8, y = 0.1, just = c("left", "top"), fontfamily = "Helvetica",
         fontsize = 14, fontface = "bold")
plotGG(pappa_age_boxplot, x = 2, width = 2.5, height = 1.6, y = 0)

## C - PAPPA GWAS variant eQTL boxplot
plotText("C", x = 0.1, y = 1.6, just = c("left", "top"), fontfamily = "Helvetica",
         fontsize = 14, fontface = "bold")

plotGG(pappa_boxplot_gwas, x = 0.5, y = 1.65, width = 3.8, height = 1.5)
plotText(label = paste0("risk allele: ", gwas_pappa$EA), fontsize = 6, fontfamily = "Helvetica",
         x = 2.5, y = 1.925)

## D - PAPPA QTL signal, Hi-C, signals, and gene
plotText("D", x = 0.1, y = 3.2, just = c("left", "top"), fontfamily = "Helvetica",
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
                            y = 4.4, height = 1, 
                            sigline = FALSE,just = c("left", "bottom"),
                            fill = colorby("LDgrp",
                                           palette = colorRampPalette(c("#262C74",
                                                                        "#98CDED",
                                                                        "#499A53",
                                                                        "#EEA741",
                                                                        "#DD3931",
                                                                        "grey"))),
                            snpHighlights = data.frame(snp = c(pappa_pbs_rsid, gwas_pappa$rsID),
                                                       pch = c(24, 23),
                                                       cex = c(0.75, 0.75),
                                                       col = c("black", "black")))
annoYaxis(plot = gwas_pappa_plot, at = seq(0, gwas_pappa_ylim, 4), 
          axisLine = TRUE, fontsize = 6, gp = gpar(fontfamily = "Helvetica"))

plotText(
  label = "-log10(p-value)", x = 0.275, y = 3.9, rot = 90,
  fontsize = 6, just = "center",
  default.units = "inches", fontfamily = "Helvetica"
)
plotText(label = "Total Hip Replacement\nGWAS", x = 4.45, y = 3.35, 
         just = c("right", "top"), fontface = "bold",
         fontsize = 8, fontfamily = "Helvetica", lineheight = 0.8)

# LD legend
grid.points(x = 0.7, y = 3.45, default.units = "native", pch = 23, size = unit(0.5, "char"))
plotText(label = paste0(gwas_pappa$rsID, " (GWAS index)"),
         fontsize = 7, fontfamily = "Helvetica",
         just = "left", x = 0.8, y = 3.45)
plotLegend(legend = c("0.8 - 1.0",
                      "0.6 - 0.8",
                      "0.4 - 0.6",
                      "0.2 - 0.4",
                      "0.0 - 0.2"),
           fill = c("#DD3931", "#EEA741", "#499A53","#98CDED","#262C74"),
           x = 0.75, y = 3.65, width = 0.1, height = 0.4, border = FALSE, 
           fontsize = 6)
plotText(label = "r2", x = 1.375, y = 3.6, fontsize = 6)
# PBS eQTL
pbs_pappa <- plotManhattan(data = pappa_signals |> filter(Condition == "PBS"),
                           params = pappa_region, 
                           range = c(0, qtl_pappa_ylim),
                           snpHighlights = data.frame(snp = c(gwas_pappa$rsID,
                                                              pappa_pbs_rsid),
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
                           y = 4.45, height = 1)
annoYaxis(plot = pbs_pappa, at = seq(0, qtl_pappa_ylim, 2), 
          axisLine = TRUE, fontsize = 6, gp = gpar(fontfamily = "Helvetica"))
plotText(
  label = "-log10(p-value)", x = 0.275, y = 5, rot = 90,
  fontsize = 6, just = "center",
  default.units = "inches", fontfamily = "Helvetica"
)
plotText(label = "eQTL", x = 4.45, y = 4.55, just = "right", fontface = "bold",
         fontsize = 8, fontfamily = "Helvetica")
plotText(label = paste0("PP4 = ", 
                        signif(pappa_colocs |> filter(Condition == "PBS") |> pull(PP4_pbs), digits = 3)), 
         x = 4.45, y = 4.65, 
         just = "right", fontsize = 7, fontfamily = "Helvetica")


# LD legend
grid.points(x = 0.7, y = 4.635, default.units = "native", pch = 24, size = unit(0.45, "char"))
plotText(label = paste0(pappa_pbs_rsid, " (eQTL index)"),
         fontsize = 7, fontfamily = "Helvetica",
         just = "left", x = 0.8, y = 4.625)
plotLegend(legend = c("0.8 - 1.0",
                      "0.6 - 0.8",
                      "0.4 - 0.6",
                      "0.2 - 0.4",
                      "0.0 - 0.2"),
           fill = c("#DD3931", "#EEA741", "#499A53","#98CDED","#262C74"),
           x = 0.75, y = 4.8, width = 0.1, height = 0.4, border = FALSE, 
           fontsize = 6)
plotText(label = "r2", x = 1.375, y = 4.75, fontsize = 6)

# Hi-C (static loop in PBS and FN-f)
hic_plot <- plotHicRectangle(data = "data/hic/CQTL_inter_30.hic",
                 resolution = 10000, params = pappa_region,
                 zrange = c(0, 130),
                 height = 1.15, y = 5.5)
annoHeatmapLegend(plot = hic_plot, x = 0.5, y = 6.075, width = 0.075, height = 0.7,
                  fontsize = 6, fontcolor = "black", just = "right")
annoPixels(plot = hic_plot, data = pappa_loop_10kb, params = pappa_region, shift = 5,
           type = "arrow")

# ATAC
plotText(label = "ATAC", rot = 90, fontcolor = yl_gn_bu[4],
         fontsize = 7, x = 0.45, y = 6.85, fontfamily = "Helvetica")
plotMultiSignal(data = list("data/atac/CQTL_ATAC_AM7754_AM7755_AM7763_PBS.bw",
                            "data/atac/CQTL_ATAC_AM7754_AM7755_AM7763_FNF.bw"),
           params = pappa_region, 
           y = 6.85, height = 0.35, linecolor = yl_gn_bu[4], just = c("left", "bottom"),
           fill = yl_gn_bu[4], gapdistance = 0.05)
plotText(label = "PBS", fontcolor = yl_gn_bu[4],
         fontsize = 5, x = 0.55, y = 6.7, just = "left", fontfamily = "Helvetica")
plotText(label = "FN-f", fontcolor = yl_gn_bu[4],
         fontsize = 5, x = 0.55, y = 6.9, just = "left", fontfamily = "Helvetica")

# CNR
plotText(label = "H3K27ac", rot = 90, fontcolor = yl_gn_bu[6],
         fontsize = 7, x = 0.45, y = 7.25, fontfamily = "Helvetica")
plotMultiSignal(data = list("data/cnr/CHON_CHON_PBS.bw",
                            "data/cnr/CHON_CHON_FNF.bw"),
                params = pappa_region, 
                y = 7.25, height = 0.35, linecolor = yl_gn_bu[6], just = c("left", "bottom"), 
                fill = yl_gn_bu[6], gapdistance = 0.05)
plotText(label = "PBS", fontcolor = yl_gn_bu[6],
         fontsize = 5, x = 0.55, y = 7.1, just = "left", fontfamily = "Helvetica")
plotText(label = "FN-f", fontcolor = yl_gn_bu[6],
         fontsize = 5, x = 0.55, y = 7.29, just = "left", fontfamily = "Helvetica")

# RNA 
plotText(label = "RNA", rot = 90, fontcolor = yl_gn_bu[8],
         fontsize = 7, x = 0.45, y = 7.7, fontfamily = "Helvetica")
plotMultiSignal(data = list("data/rna/CQTL_CTL.bw",
                            "data/rna/CQTL_FNF.bw"),
                params = pappa_region, 
                y = 7.65, height = 0.35, linecolor = yl_gn_bu[8], 
                fill = yl_gn_bu[8], gapdistance = 0.05, just = c("left", "bottom"))
plotText(label = "PBS", fontcolor = yl_gn_bu[8],
         fontsize = 5, x = 0.55, y = 7.5, just = "left", fontfamily = "Helvetica")
plotText(label = "FN-f", fontcolor = yl_gn_bu[8],
         fontsize = 5, x = 0.55, y = 7.7, just = "left", fontfamily = "Helvetica")
# Genes
pappa_gene_plot <- plotGenes(params = pappa_region, 
                             geneHighlights = data.frame("gene" = "PAPPA",
                                                         "color" = yl_gn_bu[6]),
          y = 7.85, height = 0.5)
annoGenomeLabel(plot = pappa_gene_plot, y = 8.35,
                x = 0.55, fontsize = 8)

# Loop anchors
annoHighlight(plot = pappa_gene_plot, params = pappa_region,
              chromstart = start1(pappa_loop_30kb),
              chromend = end1(pappa_loop_30kb),
              fill = "grey80", y = 5.5, height = 2.85)
annoHighlight(plot = pappa_gene_plot, params = pappa_region,
              chromstart = start2(pappa_loop_30kb),
              chromend = end2(pappa_loop_30kb),
              fill = "grey80", y = 3.4, height = 4.95)

dev.off()
