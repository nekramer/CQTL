library(DESeq2)
library(ComplexHeatmap)
library(tidyverse)
library(circlize)
library(plotgardener)
source("scripts/plotting_utils.R")

# Union of sig untreated and treated genes --------------------------------

ctl_sig_genes <- read_csv("data/ctl_sexDE_pval01.csv")
fnf_sig_genes <- read_csv("data/fnf_sexDE_pval01.csv")

union_sig_genes <- union(ctl_sig_genes$gene_id, fnf_sig_genes$gene_id)

# Untreated ---------------------------------------------------------------
load("data/dds_sex_ctl.rda")

# Normalized counts
dds_sex_ctl_norm <- vst(dds_sex_ctl)

## Read in data from ctl selecting union sig genes and order by log2FC
untreated_sex_degenes <- read_csv("data/ctl_sex_shrink.csv") %>%
  filter(gene_id %in% union_sig_genes) %>%
  mutate(log2FC_dir = ifelse(log2FoldChange < 0, "-", "+")) %>%
  mutate(log2FC_dir = factor(log2FC_dir, levels = c("-", "+"))) %>%
  # Split - and + 
  arrange(log2FC_dir) %>%
  # Arrange within each log2FC_dir group
  group_by(log2FC_dir) %>%
  arrange(desc(abs(log2FoldChange)), .by_group = TRUE) %>%
  ungroup()

# Subset norm counts for significant untreated genes
sexnormCounts_untreated <- assay(dds_sex_ctl_norm[untreated_sex_degenes$gene_id,]) %>% 
  as.data.frame()
sexnormCounts_untreated <- sexnormCounts_untreated[match(untreated_sex_degenes$gene_id, rownames(sexnormCounts_untreated)),]

# Reorder into M and F
sexnormCounts_untreated <- sexnormCounts_untreated %>%
  dplyr::select(ends_with(c("_M", "_F"))) 

# Scale counts
untreated_sexmat_scaled <- t(apply(sexnormCounts_untreated, 1, scale))
colnames(untreated_sexmat_scaled) <- colnames(sexnormCounts_untreated)

# Age, Sex, and Race Clusters
annotations <- as.data.frame(colData(dds_sex_ctl)[,c("Age_group", "Race", "Sex")])
# Put in same order as matrix
annotations <- annotations[match(colnames(untreated_sexmat_scaled), rownames(annotations)),]

annotationObjects <- HeatmapAnnotation(
  df = annotations,
  col = list(Race = c("ARAB" = raceColors[1], 
                      "ASIAN" = raceColors[2],
                      "BL" = raceColors[3], 
                      "C" = raceColors[4],
                      "HISP" = raceColors[5],
                      "Unknown" = "grey"),
             Age_group = c("30-39" = ageColors[1],
                     "40-49" = ageColors[2],
                     "50-59" = ageColors[3],
                     "60-69" = ageColors[4],
                     "70-79" = ageColors[5],
                     "80-89" = ageColors[6]),
             Sex = sexColors),
  annotation_name_gp = gpar(fontfamily = "Helvetica",
                            fontsize = 8),
  annotation_label = c("Age", "Race", "Sex"),
  which = "column"
  )

cluster_annotations <- untreated_sex_degenes %>% 
  dplyr::select(gene_id, log2FC_dir) %>%
  column_to_rownames(var = "gene_id")

clusters <- HeatmapAnnotation(
  df = cluster_annotations,
  col = list(log2FC_dir = c("-" = sexColors[["F"]], 
                         "+" = sexColors[["M"]])),
  annotation_name_gp = gpar(fontfamily = "Helvetica",
                            fontface = "bold",
                            fontsize = 0),
  which = "row",
  simple_anno_size = unit(3, "mm")
)

sex_heatmap_control <- pheatmap(untreated_sexmat_scaled, 
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



sex_heatmapLegend <- Legend(at = c(-3, 3),
                        col_fun = colorRamp2(breaks = seq(-3, 3),
                                             colors = heatmapColors),
                        border = NA,
                        title_gp = gpar(fontsize = 0),
                        labels_gp = gpar(fontfamily = "Helvetica",
                                         fontsize = 8),
                        legend_width = unit(4.55, "in"),
                        grid_height = unit(0.11, "in"),
                        direction = "horizontal"
)

sex_heatmap_controlGrob <- grid.grabExpr(draw(sex_heatmap_control,
                                  show_annotation_legend = FALSE,
                                  show_heatmap_legend = FALSE,
                                  background = "transparent"))
sex_heatmapLegendGrob <- grid.grabExpr(draw(sex_heatmapLegend))

# Plot layout together with plotgardener ----------------------------------

pdf(file = "plots/DE_heatmap_sex.pdf", width = 5.75, height = 5.25)
pageCreate(width = 5.75, height = 5.25, showGuides = FALSE)

plotGG(plot = sex_heatmap_controlGrob, x = 0, y = 0, 
       height = 5, width = 5)


# Colorbar
plotGG(plot = sex_heatmapLegendGrob,  x = 0.075, y = 5,
       width = 4.55, height = 0.11)
# Colorbar title
plotText(label = "Relative Expression", fontfamily = "Helvetica",
         fontsize = 8, x = 2.375, y = 5.1, just = "top")


# Age legend
plotRect(x = 5, y = 0.15, width = unit(3, "mm"), 
         height = unit(1, "mm"), linecolor = NA, fill = ageColors[6],
         just = "left")
plotRect(x = unit(5, "in") + unit(3, "mm"), y = 0.15, width = unit(3, "mm"), 
         height = unit(1, "mm"), linecolor = NA, fill = ageColors[5],
         just = "left")
plotRect(x = unit(5, "in") + unit(3*2, "mm"), y = 0.15, width = unit(3, "mm"), 
         height = unit(1, "mm"), linecolor = NA, fill = ageColors[4],
         just = "left")
plotRect(x = unit(5, "in") + unit(3*3, "mm"), y = 0.15, width = unit(3, "mm"), 
         height = unit(1, "mm"), linecolor = NA, fill = ageColors[3],
         just = "left")
plotRect(x = unit(5, "in") + unit(3*4, "mm"), y = 0.15, width = unit(3, "mm"), 
         height = unit(1, "mm"), linecolor = NA, fill = ageColors[2],
         just = "left")
plotRect(x = unit(5, "in") + unit(3*5, "mm"), y = 0.15, width = unit(3, "mm"), 
         height = unit(1, "mm"), linecolor = NA, fill = ageColors[1],
         just = "left")
plotText(label = "30", 
         x = 5,
         y = 0.18,
         fontfamily = "Helvetica", fontsize = 4, just = c("left", "top"))
plotText(label = "40", 
         x = unit(5, "in") + unit(3, "mm"),
         y = 0.18,
         fontfamily = "Helvetica", fontsize = 4, just = c("left", "top"))
plotText(label = "50", 
         x = unit(5, "in") + unit(3*2, "mm"),
         y = 0.18,
         fontfamily = "Helvetica", fontsize = 4, just = c("left", "top"))
plotText(label = "60", 
         x = unit(5, "in") + unit(3*3, "mm"),
         y = 0.18,
         fontfamily = "Helvetica", fontsize = 4, just = c("left", "top"))
plotText(label = "70", 
         x = unit(5, "in") + unit(3*4, "mm"),
         y = 0.18,
         fontfamily = "Helvetica", fontsize = 4, just = c("left", "top"))
plotText(label = "80", 
         x = unit(5, "in") + unit(3*5, "mm"),
         y = 0.18,
         fontfamily = "Helvetica", fontsize = 4, just = c("left", "top"))

# Race Legend
plotRect(x = unit(5, "in") + unit(3*6, "mm"), 
         y = 0.35, width = unit(3, "mm"), 
         height = unit(1, "mm"), linecolor = NA, fill = raceColors[5],
         just = "right")
plotRect(x = unit(5, "in") + unit(3*5, "mm"), 
         y = 0.35, width = unit(3, "mm"), 
         height = unit(1, "mm"), linecolor = NA, fill = raceColors[4],
         just = "right")
plotRect(x = unit(5, "in") + unit(3*4, "mm"), 
         y = 0.35, width = unit(3, "mm"), 
         height = unit(1, "mm"), linecolor = NA, fill = raceColors[3],
         just = "right")
plotRect(x = unit(5, "in") + unit(3*3, "mm"), 
         y = 0.35, width = unit(3, "mm"), 
         height = unit(1, "mm"), linecolor = NA, fill = raceColors[2],
         just = "right")
plotRect(x = unit(5, "in") + unit(3*2, "mm"), 
         y = 0.35, width = unit(3, "mm"), 
         height = unit(1, "mm"), linecolor = NA, fill = raceColors[1],
         just = "right")
plotText(label = "HI", x = unit(5, "in") + unit(3*5.5, "mm"),
         y = 0.425,
         fontsize = 5, fontfamily = "Helvetica")
plotText(label = "C", x = unit(5, "in") + unit(3*4.5, "mm"),
         y = 0.425,
         fontsize = 5, fontfamily = "Helvetica")
plotText(label = "BL", x= unit(5, "in") + unit(3*3.5, "mm"),
         y = 0.425,
         fontsize = 5, fontfamily = "Helvetica")
plotText(label = "AS", x= unit(5, "in") + unit(3*2.5, "mm"),
         y = 0.425,
         fontsize = 5, fontfamily = "Helvetica")
plotText(label = "AR", x= unit(5, "in") + unit(3*1.5, "mm"),
         y = 0.425,
         fontsize = 5, fontfamily = "Helvetica")


# Sex legend
plotRect(x = unit(5, "in") + unit(3*6, "mm"), 
         y = 0.59, width = unit(3, "mm"), 
         height = unit(1, "mm"), linecolor = NA, fill = sexColors[["M"]],
         just = "right")
plotRect(x = unit(5, "in") + unit(3*5, "mm"), 
         y = .59, width = unit(3, "mm"), 
         height = unit(1, "mm"), linecolor = NA, fill = sexColors[["F"]],
         just = "right")
plotText(label = "M", x = unit(5, "in") + unit(3*5.5, "mm"),
         y = 0.655,
         fontsize = 5, fontfamily = "Helvetica")
plotText(label = "F", x = unit(5, "in") + unit(3*4.5, "mm"),
         y = 0.655,
         fontsize = 5, fontfamily = "Helvetica")

dev.off()