library(DESeq2)
library(ComplexHeatmap)
library(tidyverse)
library(circlize)
library(patchwork)
library(ggtext)
library(ggstar)
library(plotgardener)
library(googledrive)
library(googlesheets4)
source("../utils.R")
source("../plotting_utils.R")

# SEX DE HEATMAP -----------------------------------------------------------------

ctl_sig_genes <- read_csv("data/sex_de/ctl_sexDE_pval01.csv")
fnf_sig_genes <- read_csv("data/sex_de/fnf_sexDE_pval01.csv")

union_sig_genes <- union(ctl_sig_genes$gene_id, fnf_sig_genes$gene_id)

load("data/sex_de/dds_sex_ctl.rda")

# Normalized counts
dds_sex_ctl_norm <- vst(dds_sex_ctl)

## Read in data from ctl selecting union sig genes and order by log2FC
untreated_sex_degenes <- read_csv("data/sex_de/ctl_sex_shrink.csv") |> 
  filter(gene_id %in% union_sig_genes) |> 
  mutate(log2FC_dir = ifelse(log2FoldChange < 0, "-", "+")) |> 
  mutate(log2FC_dir = factor(log2FC_dir, levels = c("-", "+"))) |> 
  # Split - and + 
  arrange(log2FC_dir) |> 
  # Arrange within each log2FC_dir group
  group_by(log2FC_dir) |> 
  arrange(desc(abs(log2FoldChange)), .by_group = TRUE) |> 
  ungroup()

# Subset norm counts for significant untreated genes
sexnormCounts_untreated <- assay(dds_sex_ctl_norm[untreated_sex_degenes$gene_id,]) |>  
  as.data.frame()
sexnormCounts_untreated <- sexnormCounts_untreated[match(untreated_sex_degenes$gene_id, rownames(sexnormCounts_untreated)),]

# Reorder into M and F
sexnormCounts_untreated <- sexnormCounts_untreated |> 
  dplyr::select(ends_with(c("_M", "_F"))) 

# Scale counts
untreated_sexmat_scaled <- t(apply(sexnormCounts_untreated, 1, scale))
colnames(untreated_sexmat_scaled) <- colnames(sexnormCounts_untreated)

# Age, Sex, and Race Clusters
annotations <- as.data.frame(colData(dds_sex_ctl)[,c("Age_group", "Ancestry", "Sex")])
# Put in same order as matrix
annotations <- annotations[match(colnames(untreated_sexmat_scaled), rownames(annotations)),]

annotationObjects <- ComplexHeatmap::HeatmapAnnotation(
  df = annotations,
  col = list(Ancestry = c("AFR" = ancestryColors[1], 
                      "AMR" = ancestryColors[2],
                      "EAS" = ancestryColors[3], 
                      "EUR" = ancestryColors[4],
                      "SAS" = ancestryColors[5]),
             Age_group = c("31-40" = ageColors[1],
                           "41-50" = ageColors[2],
                           "51-60" = ageColors[3],
                           "61-70" = ageColors[4],
                           "71-80" = ageColors[5],
                           "81-90" = ageColors[6]),
             Sex = sexColors),
  annotation_name_gp = gpar(fontfamily = "Helvetica",
                            fontsize = 6),
  annotation_label = c("Age", "Ancestry", "Sex"),
  which = "column",
  simple_anno_size = unit(3, 'mm')
)

cluster_annotations <- untreated_sex_degenes |> 
  dplyr::select(gene_id, log2FC_dir) |> 
  column_to_rownames(var = "gene_id")

clusters <- ComplexHeatmap::HeatmapAnnotation(
  df = cluster_annotations,
  col = list(log2FC_dir = c("-" = sexColors[["F"]], 
                            "+" = sexColors[["M"]])),
  annotation_name_gp = gpar(fontfamily = "Helvetica",
                            fontface = "bold",
                            fontsize = 0),
  which = "row",
  simple_anno_size = unit(3, "mm")
)

sex_heatmap_control <- ComplexHeatmap::pheatmap(untreated_sexmat_scaled, 
                                show_rownames = FALSE,
                                top_annotation = annotationObjects,
                                right_annotation = clusters,
                                row_title = NULL,
                                cluster_cols = FALSE,
                                cluster_rows = FALSE,
                                show_colnames = FALSE,
                                color = heatmapColors,
                                breaks = seq(-3, 3),
                                show_row_dend = FALSE,
                                show_column_dend = FALSE)

sex_heatmapLegend <- ComplexHeatmap::Legend(at = c(-3, 3),
                            col_fun = colorRamp2(breaks = seq(-3, 3),
                                                 colors = heatmapColors),
                            border = NA,
                            title_gp = gpar(fontsize = 0),
                            labels_gp = gpar(fontfamily = "Helvetica",
                                             fontsize = 6),
                            legend_width = unit(2.48, "in"),
                            grid_height = unit(0.09, "in"),
                            direction = "horizontal"
)

sex_heatmap_controlGrob <- grid.grabExpr(draw(sex_heatmap_control,
                                              show_annotation_legend = FALSE,
                                              show_heatmap_legend = FALSE,
                                              background = "transparent"))
sex_heatmapLegendGrob <- grid.grabExpr(draw(sex_heatmapLegend))

save(sex_heatmap_controlGrob, file = "plots/sexageDE_Fig2/sex_heatmap_controlGrob.rda")
save(sex_heatmapLegendGrob, file = "plots/sexageDE_Fig2/sex_heatmapLegendGrob.rda")

# SEX DE OVERLAP WITH SEX-BIASED GENES FROM GTEX ---------------------------------

## Wrangle data for heatmap
# Significant from PBS and FNF
ctl_sig_genes <- read_csv("data/sex_de/ctl_sexDE_pval01.csv", 
                          col_select = c("gene_id", "symbol", 
                                         "log2FoldChange", "lfcSE"))
fnf_sig_genes <- read_csv("data/sex_de/fnf_sexDE_pval01.csv",
                          col_select = c("gene_id", "symbol", 
                                         "log2FoldChange", "lfcSE"))
union_sig_genes <- union(ctl_sig_genes$gene_id, fnf_sig_genes$gene_id)

sex_genes <- bind_rows(ctl_sig_genes |> 
                         filter(gene_id %in% union_sig_genes),
                       fnf_sig_genes |> 
                         filter(gene_id %in% union_sig_genes)) |> 
  distinct(gene_id, .keep_all = TRUE) |> 
  mutate(sex = ifelse(log2FoldChange < 0, "female", "male"))

# Pull union list of our data
chond_sbgenes <- sex_genes |> 
  mutate(sex = ifelse(log2FoldChange < 0, "female", "male")) |> 
  mutate(tissue = "Chondrocytes") |> 
  dplyr::rename(effsize = log2FoldChange) |> 
  dplyr::rename(effsize_se = lfcSE)

# Read in gtex sex-biased genes 
gtex_signif_sbgenes <- 
  read_delim("data/GTEx/GTEx_Analysis_v8_sbgenes/signif.sbgenes.txt") |> 
  # convert IDs to ones compatible with ours
  mutate(gene_id = gsub("\\..*", "", gene)) |> 
  # Flip sign of effsize to match ours (gtex positive is female and negative is male)
  mutate(effsize = -1*effsize) |> 
  mutate(sex = ifelse(effsize < 0, "female", "male")) |> 
  # filter for genes in chond_sbgenes
  filter(gene_id %in% chond_sbgenes$gene_id) |> 
  mutate(tissue = gsub("_", " ", tissue))

# Join with our data
all_data_sbgenes <- bind_rows(chond_sbgenes |> 
                                dplyr::select(gene_id, tissue, sex, effsize, effsize_se),
                              gtex_signif_sbgenes |>  
                                dplyr::select(gene_id, tissue, sex, effsize, effsize_se)) |> 
  complete(gene_id, tissue) |> 
  left_join(chond_sbgenes |> dplyr::select(gene_id, symbol), by = "gene_id") 

# Get order of genes by number of overlaps between tissues
geneOverlaps <- all_data_sbgenes |> 
  group_by(symbol) |> 
  filter(!is.na(effsize)) |> 
  summarize(nOverlap = dplyr::n()) |> 
  arrange(nOverlap)

# Join gene ordering with merged gtex/pbs/fnf sb gene data
all_data_sbgenes <- left_join(all_data_sbgenes, geneOverlaps, by = "symbol") |> 
  group_by(nOverlap) |> 
  # Make all nOverlap groups consistent with male/female ordering
  arrange(sex, .by_group = TRUE) |> 
  ungroup()
 
# Set factors so chondrocytes are in the first column and gene symbols are ordered
# by number of overlaps
all_data_sbgenes$tissue = factor(all_data_sbgenes$tissue, levels = c("Chondrocytes", 
                                                    unique(gtex_signif_sbgenes$tissue)))
all_data_sbgenes$symbol <- factor(all_data_sbgenes$symbol, levels = unique(all_data_sbgenes$symbol))
                                 
all_data_sbgenes <- all_data_sbgenes |> 
  mutate(symbol_order_temp = factor(symbol, levels = rev(levels(symbol)))) |> 
  # Assign gene symbol and tissues numbers based on factor ordering
  mutate(tissue_no = as.numeric(tissue),
         gene_no = as.numeric(symbol_order_temp)) |> 
  dplyr::select(-symbol_order_temp)

## Plot heatmap of M vs F for each gene across tissues
sex_cell_heatmap <- ggplot(all_data_sbgenes, aes(x = tissue, y = symbol)) +
  scale_x_discrete(position = "top") +
  coord_cartesian(clip = "off") +
  geom_tile(aes(fill = sex), color = "grey20", linewidth = 0.005) +
  geom_rect(data = tibble(xmin = 0.5, xmax = 1.5, ymin = 0.5,
                          ymax = length(levels(all_data_sbgenes$symbol))+0.5),
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = NA, color = "black", inherit.aes = FALSE, linewidth = 0.25) +
  geom_rect(data = tibble(xmin = 0.5, xmax = length(levels(all_data_sbgenes$tissue)) + 0.5,
                          ymin = 0.5, 
                          ymax = (max(all_data_sbgenes$gene_no)+1 - all_data_sbgenes |> 
                            filter(nOverlap == 1) |> 
                            pull(gene_no) |> 
                            min()) + 0.5),
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = NA, color = "black", inherit.aes = FALSE, linewidth = 0.25) +
  scale_fill_manual(values = c(sexColors[["F"]], sexColors[["M"]]), na.value = NA) +
  scale_y_discrete(expand = c(0, 0)) +
  theme(panel.background = element_rect(fill = "transparent", color = "transparent"),
        plot.background = element_rect(fill = "transparent", color = "transparent"),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title = element_blank(),
        text = element_text(family = "Helvetica"),
        axis.ticks.x = element_line(linewidth = 0.25),
        axis.ticks.y = element_line(linewidth = 0.2),
        axis.ticks.length = unit(1.5, "mm"),
        legend.position = "none")
save(sex_cell_heatmap, file = "plots/sexageDE_Fig2/sex_cell_heatmap.rda")

## Write genes and tissues to table
# csv
all_data_sbgenes |> 
  arrange(tissue_no, gene_no) |> 
  dplyr::select(tissue_no, tissue, gene_no, symbol, gene_id, sex, effsize, effsize_se) |> 
  write_csv("tables/SupTable4.csv") 
# google drive
ss <- gs4_create(name = "SupTable4")

write_sheet(all_data_sbgenes |>
              arrange(tissue_no, gene_no) |> 
              dplyr::select(tissue_no, tissue, gene_no, symbol, gene_id, sex, effsize, effsize_se),
            ss, sheet = "Sheet1")

drive_mv(file = "SupTable4", path = as_dribble("CQTL paper/Figures and Tables"),
         overwrite = TRUE)

## Wrangle data for sex log2FC of chondrocyte-specific sb genes vs non chond specific 

# Get log2FC for all sig sex genes in PBS and FNF
ctl_sex_results <- read_csv("data/sex_de/ctl_sex_shrink.csv", 
                            col_select = c("gene_id", "symbol", "log2FoldChange", "lfcSE")) |> 
  filter(gene_id %in% sex_genes$gene_id) |> 
  mutate(condition = "PBS") |> 
  mutate(sex = ifelse(log2FoldChange < 0, "female", "male"))
fnf_sex_results <- read_csv("data/sex_de/fnf_sex_shrink.csv",
                            col_select = c("gene_id", "symbol", "log2FoldChange", "lfcSE")) |> 
  filter(gene_id %in% sex_genes$gene_id) |> 
  mutate(condition = "FN-f") |> 
  mutate(sex = ifelse(log2FoldChange < 0, "female", "male"))

pbs_fnf_sex_l2fc <- bind_rows(ctl_sex_results,
                             fnf_sex_results) |> 
  # Mark as chond specific or not
  mutate(group = ifelse(gene_id %in% gtex_signif_sbgenes$gene_id, "shared", "chondrocyte-<br>specific")) |> 
  mutate(condition = factor(condition, levels = c("PBS", "FN-f")),
         sex = factor(sex, levels = c("female", "male")))

## PBS/FNF chond effect sizes
chond_sb_l2fc <- ggplot(pbs_fnf_sex_l2fc, aes(x = group, y = log2FoldChange, fill = sex)) +
  geom_hline(yintercept = 0, lty = 2, linewidth = 0.2) +
  facet_wrap(vars(condition), scales = "free", nrow = 2, strip.position = "left") +
  geom_boxplot(outlier.shape = NA, color = "grey20", linewidth = 0.1) +
  scale_fill_manual(values = c(sexColors[["F"]], sexColors[["M"]])) +
  scale_y_continuous(name = "log~2~Fold Change",
                     limits = c(-2, 12), breaks = c(-2, seq(0, 12, 4)),
                     expand = c(0,0)) +
  theme(axis.title.y = element_markdown(family = "Helvetica", size = 5,
                                        margin = margin(r = -1)),
        axis.title.x = element_blank(),
        legend.position = "none",
        panel.background = element_rect(fill = "transparent", color = "transparent"),
        plot.background = element_rect(fill = "transparent", color = "transparent"),
        strip.background = element_blank(),
        strip.text.y.left = element_text(family = "Helvetica", size = 5, margin = margin(r=0)),
        strip.placement = "outside",
        axis.ticks.x = element_blank(),
        axis.line.y = element_line(color = "black", 
                                   linewidth = 0.2),
        axis.ticks.length.y = unit(-0.1, "cm"),
        axis.ticks.y = element_line(color = "black",
                                    linewidth = 0.2),
        axis.text.y = element_markdown(family = "Helvetica", size = 5, color = "black"),
        axis.text.x = element_markdown(family = "Helvetica", size = 4.75, color = "black"),
        axis.line.x = element_line(color = "grey25", linewidth = 0.25),
        axis.ticks = element_blank(),
        text = element_text(family = "Helvetica"))

save(chond_sb_l2fc, file = "plots/sexageDE_Fig2/chond_sb_l2fc.rda")

# AGE DE HEATMAP ----------------------------------------------------------

# Union of untreated and treated genes
ctl_cluster_pval05 <- read_csv("data/age_de/ctl_age_pval05clusters.csv") |> 
  distinct(gene_id, .keep_all = TRUE) |> 
  dplyr::select(gene_id, cluster) 

fnf_cluster_pval05 <- read_csv("data/age_de/fnf_age_pval05clusters.csv") |> 
  distinct(gene_id, .keep_all = TRUE) |> 
  dplyr::select(gene_id, cluster)

shared_sig_genes <- ctl_cluster_pval05 |> 
  filter(gene_id %in% intersect(ctl_cluster_pval05$gene_id, 
                                fnf_cluster_pval05$gene_id))

ctl_unique_sig_genes <- ctl_cluster_pval05 |> 
  filter(!gene_id %in% shared_sig_genes$gene_id)

fnf_unique_sig_genes <- fnf_cluster_pval05 |> 
  filter(!gene_id %in% shared_sig_genes$gene_id)

# Join the genes and order by cluster
union_sig_genes <- bind_rows(shared_sig_genes, ctl_unique_sig_genes, fnf_unique_sig_genes) |> 
  mutate(cluster = factor(cluster, levels = c("-", "+"))) |> 
  arrange(desc(cluster))

ages <- read_csv("data/age_de/fnf_age_pval05clusters.csv") %>%
  pull(Age) %>%
  unique() %>%
  sort()

# Normalized counts of PBS and FNF
load("data/age_de/dds_age_ctl_lrt.rda")
dds_age_ctl_lrt_norm <- vst(dds_age_ctl_lrt)
load("data/age_de/dds_age_fnf_lrt.rda")
dds_age_fnf_lrt_norm <- vst(dds_age_fnf_lrt)

# Subset CTL norm counts for CTL significant and shared genes
ctl_cluster_counts <- assay(dds_age_ctl_lrt[c(ctl_unique_sig_genes$gene_id, shared_sig_genes$gene_id)],) |> 
  as.data.frame() 
colnames(ctl_cluster_counts) <- paste0(str_extract(colnames(ctl_cluster_counts), "AM\\d{4}"), "_",
                                       str_extract(colnames(ctl_cluster_counts), "\\d+$"))
# Subset FNF norm counts for FNF significant
fnf_cluster_counts <- assay(dds_age_fnf_lrt[c(fnf_unique_sig_genes$gene_id)],) |> 
  as.data.frame()
colnames(fnf_cluster_counts) <- paste0(str_extract(colnames(fnf_cluster_counts), "AM\\d{4}"), "_",
                                       str_extract(colnames(fnf_cluster_counts), "\\d+$"))

all_cluster_counts <- bind_rows(ctl_cluster_counts, fnf_cluster_counts)
all_cluster_counts <- all_cluster_counts[match(union_sig_genes$gene_id,
                                               rownames(all_cluster_counts)),]

# Order columns by increasing age
all_cluster_counts <- all_cluster_counts |> 
  dplyr::select(ends_with(as.character(ages)))

# Scale counts
all_cluster_counts_scaled <- t(apply(all_cluster_counts, 1, scale))
colnames(all_cluster_counts_scaled) <- colnames(all_cluster_counts)

# Age, Sex, and Ancestry Clusters
annotations <- as.data.frame(colData(dds_age_ctl_lrt)[,c("Ancestry", "Sex", "Age")]) |> 
  rownames_to_column(var = "Sample") |> 
  separate_wider_delim(col = "Sample", delim = "_", names = c(NA, "Donor", NA, NA, NA, NA, NA)) |> 
  mutate(Donor_Age = paste0(Donor, "_", Age)) |> 
  column_to_rownames(var = "Donor_Age") |> 
  dplyr::select(-Donor)

# Put in same order as matrix
annotations <- annotations[match(colnames(all_cluster_counts_scaled), rownames(annotations)),]

age_seqColors <- colorRampPalette(rev(sequential_hcl(n = 9, palette = "Mint")))(length(ages))
names(age_seqColors) <- ages

annotationObjects <- HeatmapAnnotation(
  df = annotations,
  col = list(Ancestry = c("AFR" = ancestryColors[1], 
                          "AMR" = ancestryColors[2],
                          "EAS" = ancestryColors[3], 
                          "EUR" = ancestryColors[4],
                          "SAS" = ancestryColors[5]),
             Sex = sexColors,
             Age = age_seqColors),
  annotation_name_gp = gpar(fontfamily = "Helvetica",
                            fontsize = 6),
  which = "column",
  simple_anno_size = unit(3, "mm")
)

clusters <- HeatmapAnnotation(
  df = union_sig_genes %>%
    column_to_rownames(var = "gene_id"),
  col = list(cluster = ageClusterColors),
  annotation_name_gp = gpar(fontfamily = "Helvetica",
                            fontface = "bold",
                            fontsize = 0),
  which = "row",
  simple_anno_size = unit(3, "mm")
)

age_heatmap <- ComplexHeatmap::pheatmap(all_cluster_counts_scaled,
                                        show_rownames = FALSE,
                                        top_annotation = annotationObjects,
                                        right_annotation = clusters,
                                        row_title = NULL,
                                        cluster_cols = FALSE,
                                        cluster_rows = FALSE,
                                        show_colnames = FALSE,
                                        color = heatmapColors,
                                        breaks = seq(-3, 3),
                                        show_row_dend = FALSE,
                                        show_column_dend = FALSE)

age_heatmapLegend <- Legend(at = c(-3, 3),
                            col_fun = colorRamp2(breaks = seq(-3, 3),
                                                 heatmapColors),
                            border = NA,
                            title_gp = gpar(fontsize = 0),
                            labels_gp = gpar(fontfamily = "Helvetica",
                                             fontsize = 6),
                            legend_width = unit(2.48, "in"),
                            grid_height = unit(0.09, "in"),
                            direction = "horizontal"
)

age_heatmapGrob <- grid.grabExpr(draw(age_heatmap,
                                      show_annotation_legend = FALSE,
                                      show_heatmap_legend = FALSE,
                                      background = "transparent"))
age_heatmapLegendGrob <- grid.grabExpr(draw(age_heatmapLegend))

save(age_heatmapGrob, file = "plots/sexageDE_Fig2/age_heatmapGrob.rda")
save(age_heatmapLegendGrob, file = "plots/sexageDE_Fig2/age_heatmapLegendGrob.rda")

# OVERLAP WITH FNF AND OA DIFF GENES --------------------------------------

## SEX

# Union of sex degenes significant in control and FN-f
ctl_sex_degenes <- read_csv("data/sex_de/ctl_sexDE_pval01.csv")
fnf_sex_degenes <- read_csv("data/sex_de/fnf_sexDE_pval01.csv")

sex_degenes_union <- full_join(ctl_sex_degenes |> 
                                 dplyr::select("symbol", "gene_id", "log2FoldChange") |> 
                                 mutate(condition = "ctl"),
                               fnf_sex_degenes |> 
                                 dplyr::select("symbol", "gene_id", "log2FoldChange") |> 
                                 mutate(condition = "fnf")) |> 
  distinct(gene_id, .keep_all = TRUE) |> 
  dplyr::rename(sex_log2FoldChange = log2FoldChange)

# OA genes from RAAK
raak_de_genes <- read_csv("data/RAAK/RAAK_TableS3.csv", skip = 1) |> 
  dplyr::rename(symbol = GeneSYMBOL) |> 
  filter(Pval < 0.05) |> 
  mutate(log2FoldChange = log2(FC),
         raak_dir = ifelse(log2FoldChange < 0, "down", "up"),
         raak_sig = "RAAK") |> 
  dplyr::select(-FC, -log2FoldChange,-Pval)

# OA genes from Fisch et al 2018
fisch2018_de_genes <- read_csv("data/GSE114007/1-s2.0-S1063458418313876-mmc1.csv",
                               col_names = c("symbol", "log2FoldChange", "padj"),
                               col_select = c(1, 2, 3), n_max = 12475)  |>
  filter(padj < 0.05) |> 
  mutate(fisch_dir = ifelse(log2FoldChange < 0, "down", "up"),
         fisch_sig = "Fisch") |> 
  dplyr::select(-log2FoldChange, -padj)

# OA genes from Fu et al 2021
fu2021_de_genes <- read_csv("data/GSE168505/GSE168505_deseq_res.csv", 
                            col_select = c("symbol", "log2FoldChange", "padj")) |> 
  filter(padj < 0.05) |> 
  mutate(fu_dir = ifelse(log2FoldChange < 0, "down", "up"),
         fu_sig = "Fu") |> 
  dplyr::select(-log2FoldChange, -padj)

oa_de_genes <- full_join(raak_de_genes, fisch2018_de_genes, by = "symbol") |> 
  full_join(fu2021_de_genes, by = "symbol") |> 
  rowwise() |> 
  mutate(oa_group = paste(unique(na.omit(c(raak_dir, fisch_dir, fu_dir))), collapse = "/"),
         oa_study = paste(unique(na.omit(c(raak_sig, fisch_sig, fu_sig))), collapse = "/")) |> 
  ungroup() |> 
  dplyr::select(-raak_dir, -raak_sig, -fisch_sig, -fisch_dir, -fu_sig, -fu_dir)


# All sex de genes differential in OA
sex_oa_degenes_all <- sex_degenes_union |> 
  filter(symbol %in% oa_de_genes$symbol) |> 
  left_join(oa_de_genes, by = "symbol")
write_csv(sex_oa_degenes_all, file = "data/sex_de/sex_oa_degenes_all.csv")


load("data/sex_de/dds_sex_ctl.rda")
load("data/sex_de/dds_sex_fnf.rda")

# Gene counts for sex-specific/OA genes 
all_gene_sex_counts <- list()
for (geneRow in 1:nrow(sex_oa_degenes_all)){
  gene <- sex_oa_degenes_all[geneRow, ]
  
  ctl_gene_counts <- get_gene_sex_Counts(gene = gene[["gene_id"]], dds = dds_sex_ctl) |> 
    mutate(condition = "PBS")
  fnf_gene_counts <- get_gene_sex_Counts(gene = gene[["gene_id"]], dds = dds_sex_fnf) |> 
    mutate(condition = "FN-f")
  
  all_gene_counts <- bind_rows(ctl_gene_counts, fnf_gene_counts)
  all_gene_sex_counts[[gene[["gene_id"]]]] <- all_gene_counts
}

all_sex_oa_counts <- left_join(bind_rows(all_gene_sex_counts), 
                               sex_oa_degenes_all |> 
                                 dplyr::select(symbol, gene_id, oa_group), by = "gene_id") |> 
  mutate(oa_group = factor(oa_group, levels = c("up", "down")),
         condition = factor(condition, levels = c("PBS", "FN-f")))

# Filter examples for main figure plotting
fig2_sex_oa_counts <- all_sex_oa_counts |> 
  filter(symbol %in% c("SERPINE2", "RARRES2"))


# Create data for indicating significance
gene_maxCounts <- fig2_sex_oa_counts |> 
  group_by(symbol) |> 
  summarise(maxCount = max(log2(count)))

ctl_sex_OAgenes <- read_csv("data/sex_de/ctl_sexDE_pval01.csv", 
                            col_select = c("gene_id", "symbol", "padj", "log2FoldChange")) |> 
  filter(gene_id %in% fig2_sex_oa_counts$gene_id) |> 
  mutate(effect_dir = ifelse(log2FoldChange < 0, "F", "M")) |> 
  dplyr::select(-log2FoldChange) |> 
  dplyr::rename(PBS = padj)
fnf_sex_OAgenes <- read_csv("data/sex_de/fnf_sexDE_pval01.csv", 
                            col_select = c("gene_id", "symbol", "padj", "log2FoldChange")) |> 
  filter(gene_id %in% fig2_sex_oa_counts$gene_id) |> 
  mutate(effect_dir = ifelse(log2FoldChange < 0, "F", "M")) |> 
  dplyr::select(-log2FoldChange) |> 
  dplyr::rename(`FN-f` = padj)

fig2_sex_sigline <- full_join(ctl_sex_OAgenes, fnf_sex_OAgenes) |> 
  pivot_longer(cols = c("PBS", "FN-f"), names_to = "condition", values_to = "padj") |> 
  left_join(gene_maxCounts) |> 
  filter(!is.na(padj)) |> 
  mutate(signif = "*")
fig2_sex_sigline <- merge(fig2_sex_sigline, data.frame("Sex" = c("M", "F"))) |> 
  left_join(fig2_sex_oa_counts |> dplyr::select(symbol, oa_group) |> distinct(),
            by = "symbol") |> 
  mutate(condition = factor(condition, levels = c("PBS", "FN-f")),
         oa_group = factor(oa_group, levels = c("up", "down")))


sex_oa_boxplot_examples <- ggplot(fig2_sex_oa_counts, 
                      aes(x = Sex, y = log2(count), fill = Sex)) +
  geom_boxplot(outlier.shape = NA, linewidth = 0.25) +
  ggh4x::facet_grid2(condition ~ oa_group + symbol, switch = "y", scales = "free",
                     independent = "all", labeller = as_labeller(c("up" = "Up in OA",
                                                                   "down" = "Down in OA",
                                                                   "SERPINE2" = "**SERPINE2**",
                                                                   "RARRES2" = "**RARRES2**",
                                                                   "PBS" = "PBS",
                                                                   "FN-f" = "FN-f"))) +
  geom_line(data = fig2_sex_sigline, 
            aes(x = Sex, y = maxCount + 0.5, group = symbol), inherit.aes = FALSE,
            linewidth = 0.25) +
  geom_star(data = fig2_sex_sigline,
            aes(x = 1.5, y = Inf, starshape = signif, fill = effect_dir), inherit.aes = FALSE,
            size = 1.5, starstroke = 0.25) +
  scale_fill_manual(values = c(sexColors[["F"]], sexColors[["M"]])) +
  scale_y_continuous(breaks = scales::breaks_pretty(3), 
                     name = "log~2~(normalized counts)") +
  coord_cartesian(clip = "off") +
  theme(strip.placement = "outside",
        axis.line = element_line(linewidth = 0.25),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(color = "black", linewidth = 0.25),
        axis.ticks.length.y = unit(-0.1, "cm"),
        axis.title.y = element_markdown(size = 6, family = "Helvetica",
                                        margin = margin(r = -15)),
        text = element_text(family = "Helvetica"),
        axis.text.y = element_text(color = "black", size = 6),
        axis.text.x = element_text(color = "black", size = 6, margin = margin(b = 0)),
        strip.background = element_blank(),
        strip.text.y.left = element_text(color = "black", size = 6, angle = 0,
                                         margin = margin()),
        strip.text.x.top = element_markdown(size = 8, margin = margin(t = 1)),
        panel.background = element_rect(fill = "transparent", color = "transparent"),
        plot.background = element_rect(fill = "transparent", color = "transparent"),
        panel.grid = element_blank(),
        legend.position = "none",
        panel.spacing.y = unit(0.25, "cm"))

save(sex_oa_boxplot_examples, file = "plots/sexageDE_Fig2/sex_oa_boxplot_examples.rda")

## AGE

# Get union of significant age genes
ctl_cluster_pval05 <- read_csv("data/age_de/ctl_age_pval05clusters.csv") |> 
  distinct(gene_id, .keep_all = TRUE) |> 
  dplyr::select(gene_id, symbol, cluster) |> 
  dplyr::rename(ctl_cluster = cluster)

fnf_cluster_pval05 <- read_csv("data/age_de/fnf_age_pval05clusters.csv") |> 
  distinct(gene_id, .keep_all = TRUE) |> 
  dplyr::select(gene_id, symbol, cluster) |> 
  dplyr::rename(fnf_cluster = cluster)

union_sig_genes <- full_join(ctl_cluster_pval05, fnf_cluster_pval05, 
                             by = c("gene_id", "symbol"))

# OA genes from RAAK
raak_de_genes <- read_csv("data/RAAK/RAAK_TableS3.csv", skip = 1) |> 
  dplyr::rename(symbol = GeneSYMBOL) |> 
  filter(Pval < 0.05) |> 
  mutate(log2FoldChange = log2(FC),
         raak_dir = ifelse(log2FoldChange < 0, "down", "up"),
         raak_sig = "RAAK") |> 
  dplyr::select(-FC, -log2FoldChange,-Pval)

# OA genes from Fisch et al 2018
fisch2018_de_genes <- read_csv("data/GSE114007/1-s2.0-S1063458418313876-mmc1.csv",
                               col_names = c("symbol", "log2FoldChange", "padj"),
                               col_select = c(1, 2, 3), n_max = 12475)  |>
  filter(padj < 0.05) |> 
  mutate(fisch_dir = ifelse(log2FoldChange < 0, "down", "up"),
         fisch_sig = "Fisch") |> 
  dplyr::select(-log2FoldChange, -padj)

# OA genes from Fu et al 2021
fu2021_de_genes <- read_csv("data/GSE168505/GSE168505_deseq_res.csv", 
                            col_select = c("symbol", "log2FoldChange", "padj")) |> 
  filter(padj < 0.05) |> 
  mutate(fu_dir = ifelse(log2FoldChange < 0, "down", "up"),
         fu_sig = "Fu") |> 
  dplyr::select(-log2FoldChange, -padj)

oa_de_genes <- full_join(raak_de_genes, fisch2018_de_genes, by = "symbol") |> 
  full_join(fu2021_de_genes, by = "symbol") |> 
  rowwise() |> 
  mutate(oa_group = paste(unique(na.omit(c(raak_dir, fisch_dir, fu_dir))), collapse = "/"),
         oa_study = paste(unique(na.omit(c(raak_sig, fisch_sig, fu_sig))), collapse = "/")) |> 
  ungroup() |> 
  dplyr::select(-raak_dir, -raak_sig, -fisch_sig, -fisch_dir, -fu_sig, -fu_dir)

# All age de genes differential in OA
age_oa_degenes_all <- union_sig_genes |> 
  filter(symbol %in% oa_de_genes$symbol) |> 
  left_join(oa_de_genes, by = "symbol")
write_csv(age_oa_degenes_all, file = "data/age_de/age_oa_degenes_all.csv")

# Get counts for age OA de genes in PBS and FN-f
load("data/age_de/dds_age_ctl_lrt.rda")
load("data/age_de/dds_age_fnf_lrt.rda")

all_gene_age_counts <- list()
for (geneRow in 1:nrow(age_oa_degenes_all)){
  gene <- age_oa_degenes_all[geneRow, ]
  
  ctl_gene_counts <- get_gene_age_Counts(gene = gene[["gene_id"]], 
                                         dds = dds_age_ctl_lrt) |> 
    mutate(condition = "PBS")
  fnf_gene_counts <- get_gene_age_Counts(gene = gene[["gene_id"]], 
                                         dds = dds_age_fnf_lrt) |> 
    mutate(condition = "FN-f")
  
  all_gene_counts <- bind_rows(ctl_gene_counts, fnf_gene_counts)
  all_gene_age_counts[[gene[["gene_id"]]]] <- all_gene_counts
}

# Join with OA data
all_age_oa_counts <- left_join(bind_rows(all_gene_age_counts),
                                 age_oa_degenes_all,
                                 by = "gene_id") |> 
  mutate(ctl_cluster = ifelse(is.na(ctl_cluster), fnf_cluster, ctl_cluster),
         fnf_cluster = ifelse(is.na(fnf_cluster), ctl_cluster, fnf_cluster)) |> 
  mutate(ctl_cluster = factor(ctl_cluster, levels = c("-", "+"))) |> 
  mutate(oa_group = factor(oa_group, levels = c("up", "down")),
         condition = factor(condition, levels = c("PBS", "FN-f"))) |> 
  # Split data into age groups
  mutate(Age_group = case_when(Age <= 50 ~ "<50",
                               Age >= 65 ~ "65+",
                               Age > 50 & Age < 65 ~ "50-65")) |> 
  mutate(Age_group = factor(Age_group, levels = c("<50", "50-65", "65+")))


# Filter examples for main figure plotting
fig2_age_oa_counts <- all_age_oa_counts |> 
  filter(symbol %in% c("EDA2R", "IRS1"))

# Create data for indicating significance
gene_maxCounts <- fig2_age_oa_counts |> 
  group_by(symbol) |> 
  summarise(maxCount = max(log2(count)))

fig2_age_sigline <- age_oa_degenes_all |> 
  filter(symbol %in% c("EDA2R", "IRS1")) |> 
  pivot_longer(cols = ends_with("cluster"), 
               names_to = "condition", 
               values_to = "cluster") |> 
  filter(!is.na(cluster)) |> 
  mutate(condition = gsub("_cluster", "", condition)) |> 
  mutate(condition = ifelse(condition == "ctl", "PBS", "FN-f")) |> 
  mutate(cluster = factor(cluster, levels = c("-", "+")),
         condition = factor(condition, levels = c("PBS", "FN-f")),
         oa_group = factor(oa_group, levels = c("up", "down"))) |> 
  left_join(gene_maxCounts, by = "symbol") |> 
  merge(data.frame("Age_group" = c("<50", "50-65", "65+")))


age_oa_boxplot_examples <- ggplot(fig2_age_oa_counts, aes(x = Age_group, y = log2(count))) +
  geom_boxplot(outlier.shape = NA, aes(fill = ctl_cluster),
               linewidth = 0.25)  +
  ggh4x::facet_grid2(condition ~ oa_group + symbol, switch = "y", scales = "free",
                     independent = "all", labeller = as_labeller(c("up" = "Up in OA",
                                                                   "down" = "Down in OA",
                                                                   "EDA2R" = "**EDA2R**",
                                                                   "IRS1" = "**IRS1**",
                                                                   "PBS" = "PBS",
                                                                   "FN-f" = "FN-f"))) +
  scale_fill_manual(values = c("-" = ageClusterColors[["-"]], "+" = ageClusterColors[["+"]])) +
  geom_line(data = fig2_age_sigline, 
            aes(x = Age_group, y = maxCount + 0.5, group = symbol), inherit.aes = FALSE,
            linewidth = 0.25) +
  geom_star(data = fig2_age_sigline,
            aes(x = 2, y = Inf, fill = cluster), inherit.aes = FALSE,
            size = 1.5, starstroke = 0.25) +
  scale_y_continuous(breaks = scales::breaks_pretty(3), 
                     name = "log~2~(normalized counts)") +
  coord_cartesian(clip = "off") +
  theme(strip.placement = "outside",
        axis.line = element_line(linewidth = 0.25),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(color = "black", linewidth = 0.25),
        axis.ticks.length.y = unit(-0.1, "cm"),
        axis.title.y = element_markdown(size = 6, family = "Helvetica",
                                        margin = margin(r = -15)),
        text = element_text(family = "Helvetica"),
        axis.text.y = element_text(color = "black", size = 6),
        axis.text.x = element_text(color = "black", size = 6, margin = margin(b = 0)),
        strip.background = element_blank(),
        strip.text.y.left = element_text(color = "black", size = 6, angle = 0,
                                         margin = margin()),
        strip.text.x.top = element_markdown(size = 8, margin = margin(t = 1)),
        panel.background = element_rect(fill = "transparent", color = "transparent"),
        plot.background = element_rect(fill = "transparent", color = "transparent"),
        panel.grid = element_blank(),
        legend.position = "none",
        panel.spacing.y = unit(0.5, "cm"),
        plot.title = element_text(hjust = 0.5, size = 10, margin = margin(b = -1)))

save(age_oa_boxplot_examples, file = "plots/sexageDE_Fig2/age_oa_boxplot_examples.rda")


# AGE GO TERMS ------------------------------------------------------------

age_cluster_up_go <- 
  read_delim("data/homer/homer_age_sig_cluster_up/biological_process.txt") |> 
  mutate(pval = exp(1)^logP) |> 
  filter(pval < 0.01)
page_cluster_up_go <- reduceGO(age_cluster_up_go,
                              category = "clust_up")

## Format and write to table
age_upgo_table <- age_cluster_up_go |> 
  dplyr::select(-`Entrez Gene IDs`, -pval, -logP, -size, -termUniqueness, -termUniquenessWithinCluster,
                -termDispensability, -category) |> 
  relocate(`-log10pval`, .after = Enrichment) |> 
  arrange(desc(`-log10pval`))

write_csv(age_upgo_table, file = "tables/SupTable6A.csv")

age_cluster_up_goPlotting <- age_cluster_up_go |> 
  filter(parentTerm %in% c("cartilage condensation", "humoral immune response",
                           "positive regulation of transcription by RNA polymerase II", 
                           "tryptophan catabolic process to kynurenine", "thrombin-activated receptor signaling pathway"))

age_cluster_down_go <-
  read_delim("data/homer/homer_age_sig_cluster_down/biological_process.txt") |> 
  mutate(pval = exp(1)^logP) |> 
  filter(pval < 0.01)
age_cluster_down_go <- reduceGO(age_cluster_down_go,
                                category = "clust_down") 

## Format and write to table
age_downgo_table <- age_cluster_down_go |> 
  dplyr::select(-`Entrez Gene IDs`, -pval, -logP, -size, -termUniqueness, -termUniquenessWithinCluster,
                -termDispensability, -category) |> 
  relocate(`-log10pval`, .after = Enrichment) |> 
  arrange(desc(`-log10pval`))

write_csv(age_downgo_table, file = "tables/SupTable6B.csv")

age_cluster_down_goPlotting <- age_cluster_down_go |> 
  filter(parentTerm %in% c("type I interferon-mediated signaling pathway", 
                           "regulation of response to drug", 
                           "negative regulation of signaling",
                           "negative regulation of signal transduction",
                           "negative regulation of fibroblast growth factor receptor signaling pathway"))

age_go_plotting <- bind_rows(age_cluster_up_goPlotting, age_cluster_down_goPlotting) |> 
  mutate(category = factor(category, levels = c("clust_up", "clust_down"))) |> 
  arrange(category)
age_go_plotting$parentTerm <- factor(age_go_plotting$parentTerm, 
                                     levels = rev(age_go_plotting$parentTerm))


# Plot all in barplot
age_GO_barplots <- ggplot(age_go_plotting, aes(x = `-log10pval`, y = parentTerm, fill = category)) +
  geom_bar(stat = "identity") +
  scale_x_continuous(limits = c(0, 5), expand = c(0, 0), name = "-log~10~pval",
                     breaks = seq(0, 5, 1)) +
  scale_fill_manual(values = c("clust_up" = ageClusterColors[["+"]], 
                               "clust_down" = ageClusterColors[["-"]])) +
  coord_cartesian(clip = "off") +
  facet_wrap(~category, ncol = 1, strip.position = "left", scales = "free_y") +
  geom_text(aes(x = 0, label = parentTerm), hjust = 0, family = "Helvetica",
            size = 2.75) +
  theme(panel.background = element_rect(fill = "transparent", color = "transparent"),
        plot.background = element_rect(fill = "transparent", color = "transparent"),
        text = element_text(family = "Helvetica"),
        legend.position = "None",
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_markdown(size = 6, margin = margin(t = -3)),
        axis.text.x = element_text(color = "black", size = 6,
                                   margin = margin(t = -1)),
        strip.background = element_blank(),
        strip.text.y.left = element_blank(),
        axis.ticks = element_blank(),
        axis.line.x = element_line(linewidth = 0.25),
        panel.spacing = unit(0, "mm"), 
        plot.title = element_text(hjust = 0.5, face = "bold", size = 10)) +
  ggtitle("GO Terms")

save(age_GO_barplots, file = "plots/sexageDE_Fig2/age_GO_barplots.rda")
# Assemble entire figure with plotgardener --------------------------------

pdf(file = "plots/sexageDE_Fig2/Fig2.pdf", width = 11.25, height = 6.5)

pageCreate(width = 11.25, height = 6.5, showGuides = FALSE)

## A - sex heatmap
plotText("A", x = 0.1, y = 0.2, just = c("left", "top"), fontfamily = "Helvetica",
         fontsize = 14, fontface = "bold")
plotGG(plot = sex_heatmap_controlGrob, x = 0.35, y = 0.25, 
       height = 2.75, width = 3)

# Colorbar
plotGG(plot = sex_heatmapLegendGrob, x = 0.425, y = 3,
       width = 2.48, height = 0.09)
# Colorbar title
plotText(label = "Relative Expression", fontfamily = "Helvetica",
         fontsize = 6, x = 1.65, y = 3.1, just = "top")

# Age legend
plotRect(x = 3.3, y = 0.375, width = unit(2, "mm"), 
         height = unit(1, "mm"), linecolor = NA, fill = ageColors[6],
         just = "left")
plotRect(x = unit(3.3, "in") + unit(2, "mm"), y = 0.375, width = unit(2, "mm"), 
         height = unit(1, "mm"), linecolor = NA, fill = ageColors[5],
         just = "left")
plotRect(x = unit(3.3, "in") + unit(2*2, "mm"), y = 0.375, width = unit(2, "mm"), 
         height = unit(1, "mm"), linecolor = NA, fill = ageColors[4],
         just = "left")
plotRect(x = unit(3.3, "in") + unit(2*3, "mm"), y = 0.375, width = unit(2, "mm"), 
         height = unit(1, "mm"), linecolor = NA, fill = ageColors[3],
         just = "left")
plotRect(x = unit(3.3, "in") + unit(2*4, "mm"), y = 0.375, width = unit(2, "mm"), 
         height = unit(1, "mm"), linecolor = NA, fill = ageColors[2],
         just = "left")
plotRect(x = unit(3.3, "in") + unit(2*5, "mm"), y = 0.375, width = unit(2, "mm"), 
         height = unit(1, "mm"), linecolor = NA, fill = ageColors[1],
         just = "left")
plotText(label = "31", 
         x = 3.3,
         y = 0.4,
         fontfamily = "Helvetica", fontsize = 4, just = c("left", "top"))
plotText(label = "41", 
         x = unit(3.3, "in") + unit(2, "mm"),
         y = 0.4,
         fontfamily = "Helvetica", fontsize = 4, just = c("left", "top"))
plotText(label = "51", 
         x = unit(3.3, "in") + unit(2*2, "mm"),
         y = 0.4,
         fontfamily = "Helvetica", fontsize = 4, just = c("left", "top"))
plotText(label = "61", 
         x = unit(3.3, "in") + unit(2*3, "mm"),
         y = 0.4,
         fontfamily = "Helvetica", fontsize = 4, just = c("left", "top"))
plotText(label = "71", 
         x = unit(3.3, "in") + unit(2*4, "mm"),
         y = 0.4,
         fontfamily = "Helvetica", fontsize = 4, just = c("left", "top"))
plotText(label = "81", 
         x = unit(3.3, "in") + unit(2*5, "mm"),
         y = 0.4,
         fontfamily = "Helvetica", fontsize = 4, just = c("left", "top"))

# Ancestry Legend
plotRect(x = unit(3.3, "in") + unit(3*4, "mm"), 
         y = 0.525, width = unit(3, "mm"), 
         height = unit(1, "mm"), linecolor = NA, fill = ancestryColors[5],
         just = "left")
plotRect(x = unit(3.3, "in") + unit(3*3, "mm"), 
         y = 0.525, width = unit(3, "mm"), 
         height = unit(1, "mm"), linecolor = NA, fill = ancestryColors[4],
         just = "left")
plotRect(x = unit(3.3, "in") + unit(3*2, "mm"), 
         y = 0.525, width = unit(3, "mm"), 
         height = unit(1, "mm"), linecolor = NA, fill = ancestryColors[3],
         just = "left")
plotRect(x = unit(3.3, "in") + unit(3*1, "mm"), 
         y = 0.525, width = unit(3, "mm"), 
         height = unit(1, "mm"), linecolor = NA, fill = ancestryColors[2],
         just = "left")
plotRect(x = unit(3.3, "in"), 
         y = 0.525, width = unit(3, "mm"), 
         height = unit(1, "mm"), linecolor = NA, fill = ancestryColors[1],
         just = "left")
plotText(label = "SAS", x = unit(3.3, "in") + unit(3*4.5, "mm"),
         y = 0.575,
         fontsize = 3.25, fontfamily = "Helvetica")
plotText(label = "EUR", x = unit(3.3, "in") + unit(3*3.5, "mm"),
         y = 0.575,
         fontsize = 3.25, fontfamily = "Helvetica")
plotText(label = "EAS", x = unit(3.3, "in") + unit(3*2.5, "mm"),
         y = 0.575,
         fontsize = 3.25, fontfamily = "Helvetica")
plotText(label = "AMR", x = unit(3.3, "in") + unit(3*1.5, "mm"),
         y = 0.575,
         fontsize = 3.25, fontfamily = "Helvetica")
plotText(label = "AFR", x = unit(3.3, "in") + unit(3*0.5, "mm"),
         y = 0.575,
         fontsize = 3.25, fontfamily = "Helvetica")

# Sex legend
plotRect(x = unit(3.3, "in") + unit(2, "mm"), 
         y = 0.65, width = unit(2, "mm"), 
         height = unit(1, "mm"), linecolor = NA, fill = sexColors[["M"]],
         just = "left")
plotRect(x = unit(3.3, "in"), 
         y = 0.65, width = unit(2, "mm"), 
         height = unit(1, "mm"), linecolor = NA, fill = sexColors[["F"]],
         just = "left")
plotText(label = "M", x = unit(3.3, "in") + unit(2*1.5, "mm"),
         y = 0.715,
         fontsize = 4.5, fontfamily = "Helvetica")
plotText(label = "F", x = unit(3.3, "in") + unit(2*0.5, "mm"),
         y = 0.715,
         fontsize = 4.5, fontfamily = "Helvetica")

## B - gtex-sex overlap
plotText("B", x = 4, y = 0.2, just = c("left", "top"), fontfamily = "Helvetica",
         fontsize = 14, fontface = "bold")
plotGG(plot = sex_cell_heatmap, x = 4.1, y = 0.25, width = 2.75, height = 2.9)
plotGG(plot = chond_sb_l2fc, x = 6.75, y = 1.05, width = 1.25, height = 2.15)

## C - overlap with OA, sex examples
plotText("C", x = 8, y = 0.2, just = c("left", "top"), fontfamily = "Helvetica",
         fontsize = 14, fontface = "bold")
plotGG(sex_oa_boxplot_examples, x = 8.1, y = -0.05, width = 3, height = 3.2)

## D - age heatmap
plotText("D", x = 0.1, y = 3.3, just = c("left", "top"), fontfamily = "Helvetica",
         fontsize = 14, fontface = "bold")

plotGG(plot = age_heatmapGrob, x = 0.35, y = 3.4, 
       height = 2.75, width = 3)

plotGG(plot = age_heatmapLegendGrob, x = 0.425, y = 6.15,
       width = 2.48, height = 0.09)
# Colorbar title
plotText(label = "Relative Expression", fontfamily = "Helvetica",
         fontsize = 6, x = 1.65, y = 6.225, just = "top")

# Ancestry Legend
plotRect(x = unit(3.3, "in") + unit(3*4, "mm"), 
         y = 3.5, width = unit(3, "mm"), 
         height = unit(1, "mm"), linecolor = NA, fill = ancestryColors[5],
         just = "left")
plotRect(x = unit(3.3, "in") + unit(3*3, "mm"), 
         y = 3.5, width = unit(3, "mm"), 
         height = unit(1, "mm"), linecolor = NA, fill = ancestryColors[4],
         just = "left")
plotRect(x = unit(3.3, "in") + unit(3*2, "mm"), 
         y = 3.5, width = unit(3, "mm"), 
         height = unit(1, "mm"), linecolor = NA, fill = ancestryColors[3],
         just = "left")
plotRect(x = unit(3.3, "in") + unit(3*1, "mm"), 
         y = 3.5, width = unit(3, "mm"), 
         height = unit(1, "mm"), linecolor = NA, fill = ancestryColors[2],
         just = "left")
plotRect(x = unit(3.3, "in"), 
         y = 3.5, width = unit(3, "mm"), 
         height = unit(1, "mm"), linecolor = NA, fill = ancestryColors[1],
         just = "left")
plotText(label = "SAS", x = unit(3.3, "in") + unit(3*4.5, "mm"),
         y = 3.55,
         fontsize = 3.25, fontfamily = "Helvetica")
plotText(label = "EUR", x = unit(3.3, "in") + unit(3*3.5, "mm"),
         y = 3.55,
         fontsize = 3.25, fontfamily = "Helvetica")
plotText(label = "EAS", x = unit(3.3, "in") + unit(3*2.5, "mm"),
         y = 3.55,
         fontsize = 3.25, fontfamily = "Helvetica")
plotText(label = "AMR", x = unit(3.3, "in") + unit(3*1.5, "mm"),
         y = 3.55,
         fontsize = 3.25, fontfamily = "Helvetica")
plotText(label = "AFR", x = unit(3.3, "in") + unit(3*0.5, "mm"),
         y = 3.55,
         fontsize = 3.25, fontfamily = "Helvetica")

# Sex legend
plotRect(x = unit(3.3, "in") + unit(2, "mm"), 
         y = 3.665, width = unit(2, "mm"), 
         height = unit(1, "mm"), linecolor = NA, fill = sexColors[["M"]],
         just = "left")
plotRect(x = unit(3.3, "in"), 
         y = 3.665, width = unit(2, "mm"), 
         height = unit(1, "mm"), linecolor = NA, fill = sexColors[["F"]],
         just = "left")
plotText(label = "M", x = unit(3.3, "in") + unit(2*1.5, "mm"),
         y = 3.72,
         fontsize = 4.5, fontfamily = "Helvetica")
plotText(label = "F", x = unit(3.3, "in") + unit(2*0.5, "mm"),
         y = 3.72,
         fontsize = 4.5, fontfamily = "Helvetica")

# Age legend
ages <- read_csv("data/age_de/fnf_age_pval01clusters.csv") %>%
  pull(Age) %>%
  unique() %>%
  sort()
age_seqColors <- colorRampPalette(rev(sequential_hcl(n = 9, palette = "Mint")))(length(ages))
names(age_seqColors) <- ages

xcoord <- unit(3.3, "in")
for (i in 1:length(age_seqColors)){
  plotRect(x = xcoord, y = 3.8, width = unit(18/38, "mm"),
           height = unit(1, "mm"), linecolor = NA,
           fill = age_seqColors[i], just = "left")
  xcoord <- xcoord + unit(18/38, "mm")
}

plotText(label = "34", 
         x = 3.3,
         y = 3.85,
         fontfamily = "Helvetica", fontsize = 5, just = c("left", "top"))

plotText(label = "84", 
         x = unit(3.3, "in") + unit(3*6, "mm"),
         y = 3.85,
         fontfamily = "Helvetica", fontsize = 5, just = c("right", "top"))


## E - age GO terms
plotText("E", x = 4, y = 3.3, just = c("left", "top"), fontfamily = "Helvetica",
         fontsize = 14, fontface = "bold")

plotGG(plot = age_GO_barplots, x = 4.2, y = 3.5, width = 3.75, height = 2.9)

## F - overlap with OA, sex examples
plotText("F", x = 8, y = 3.3, just = c("left", "top"), fontfamily = "Helvetica",
         fontsize = 14, fontface = "bold")

plotGG(age_oa_boxplot_examples, x = 8.1, y = 3.2, width = 3, height = 3.2)


dev.off()
