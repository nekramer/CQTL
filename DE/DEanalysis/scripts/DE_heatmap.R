library(ComplexHeatmap)
library(DESeq2)
library(tidyverse)
library(circlize)
library(RColorBrewer)
library(plotgardener)
library(grid)
source("scripts/plotting_utils.R")

load("data/differential_expression_dds.rda")
sig_degenes <- read_csv("data/sig_deGenes_pval01_l2fc2.csv") %>%
  mutate(log2FC_dir = ifelse(log2FoldChange < 0, "-", "+")) %>%
  arrange(log2FC_dir)
donorSamplesheet <- read_csv("data/donorSamplesheet.csv")

# Heatmap -----------------------------------------------------------------

# Normalized counts
dds_norm <- vst(dds)

# Grab significant genes
normCounts <- assay(dds_norm)[sig_degenes$gene_id,] %>% 
  as.data.frame() 
normCounts <- normCounts[match(sig_degenes$gene_id, rownames(normCounts)),]

# Sort columns into CTL and FNF
normCounts <- normCounts %>%
  # Sort columns into CTL and FNF
  dplyr::select(contains(c("CTL", "FNF")))

# Scale counts
mat_scaled <- t(apply(normCounts, 1, scale))
colnames(mat_scaled) <- colnames(normCounts)

# Condition, Age, Sex, Race Clusters
annotations <- as.data.frame(colData(dds)[,c("Condition", "Donor")]) %>%
  rownames_to_column(var = "Sample") %>%
  # Put in same order as matrix
  arrange(match(Sample, colnames(normCounts))) %>%
  dplyr::select(-Sample)
colnames(annotations) <- c("Condition", "Donor")
annotations <- left_join(annotations, donorSamplesheet[,c("Donor", "Sex", "Race", "Age")]) %>% 
  dplyr::select(-Donor)
annotations$Age <- cut(annotations$Age, breaks = seq(29, 91, 10),
                            labels = c("30-39", "40-49", "50-59", "60-69", "70-79", "80-89"))
annotations <- annotations %>% dplyr::select(Age, Race, Sex, Condition)

# Significant genes to highlight

downGenes <- sig_degenes %>% 
  filter(symbol %in% c("GREM1", "DKK1", "GDF5", "GDF10", "DLX5",
                       "GPX3", "COL21A1", "FZD8", "WWP2", "ALDH1A1",
                       "NOG", "MAP2K6", "SOX6", "NKX3-2", "ERG")) %>%
  arrange(log2FoldChange)

upGenes <- sig_degenes %>% 
  filter(symbol %in% c("CXCL2", "LIF", "IL6", "IL1B", "MMP13",
                       "ADAMTS4", "NFKB1", "CAMK1G", "IRAK2", "MMP1",
                       "IL17C", "CXCR4", "WNT5A", "BMP6", "COL13A1",
                       "IL11", "CRTAC1", "COL7A1", "MMP10", "CXCL1")) %>%
  arrange(desc(log2FoldChange))


annotationObjects <- HeatmapAnnotation(
  df = annotations,
  col = list(Age = c("30-39" = ageColors[6],
                     "40-49" = ageColors[5],
                     "50-59" = ageColors[4],
                     "60-69" = ageColors[3],
                     "70-79" = ageColors[2],
                     "80-89" = ageColors[1]),
             Condition = conditionColors,
             Sex = sexColors,
             Race = c("ARAB" = raceColors[1],
                      "ASIAN" = raceColors[2],
                      "BL" = raceColors[3],
                      "C" = raceColors[4],
                      "HISP" = raceColors[5])),
  annotation_name_gp = gpar(fontfamily = "Helvetica",
                            fontsize = 8),
  which = "column"
)

cluster_annotations <- sig_degenes %>% 
  dplyr::select(gene_id, log2FC_dir) %>%
  column_to_rownames(var = "gene_id")


clusters <- HeatmapAnnotation(
  df = cluster_annotations,
  col = list(log2FC_dir = log2fcColors),
  annotation_name_gp = gpar(fontfamily = "Helvetica",
                            fontface = "bold",
                            fontsize = 0),
  simple_anno_size = unit(3, "mm"),
  which = "row"
)

# h1 <- pheatmap(mat_scaled, 
#          show_rownames = FALSE,
#          top_annotation = annotationObjects,
#          right_annotation = clusters,
#          row_title =NULL,
#          cluster_cols = TRUE,
#          cluster_rows = FALSE,
#          show_colnames = FALSE,
#          color = heatmapColors,
#          breaks = seq(-3, 3),
#          show_row_dend = FALSE,
#          show_column_dend = FALSE)

# Initialize clustering of heatmap

h1 <- draw(Heatmap(mat_scaled,
              show_row_names = FALSE,
              top_annotation = annotationObjects,
              right_annotation = clusters,
              row_title = NULL,
              cluster_columns = TRUE,
              cluster_rows = FALSE,
              show_column_names = FALSE,
              col = colorRamp2(seq(-3, 3), heatmapColors),
              show_row_dend = FALSE,
              show_column_dend = FALSE))
# Get column order of clusters and swap CTL and FNF clusters
col_order <- column_order(h1)
new_col_order <- c(col_order[104:206], col_order[1:103])

# Plot heatmap with column order defined from clustering and CTL and FNF order above
h1 <- Heatmap(mat_scaled,
              show_row_names = FALSE,
              top_annotation = annotationObjects,
              right_annotation = clusters,
              row_title = NULL,
              cluster_columns = FALSE,
              cluster_rows = FALSE,
              show_column_names = FALSE,
              col = colorRamp2(seq(-3, 3), heatmapColors),
              show_row_dend = FALSE,
              show_column_dend = FALSE,
              column_order = new_col_order)

heatmapLegend <- Legend(at = c(-3, 3),
                        col_fun = colorRamp2(breaks = seq(-3, 3),
                                             colors = heatmapColors),
                        border = NA,
                        title_gp = gpar(fontsize = 0),
                        labels_gp = gpar(fontfamily = "Helvetica",
                                         fontsize = 8),
                        legend_width = unit(5.35, "in"),
                        grid_height = unit(0.11, "in"),
                        direction = "horizontal"
                        )


heatmapGrob <- grid.grabExpr(draw(h1,
                                  show_annotation_legend = FALSE,
                                  show_heatmap_legend = FALSE,
                                  background = "transparent"))
heatmapLegendGrob <- grid.grabExpr(draw(heatmapLegend))


# Plot in plotgardener layout
pdf(file = "plots/DE_heatmap.pdf", width = 6.5, height = 5.25)
pageCreate(width = 6.5, height = 5.25, showGuides = FALSE)
plotGG(plot = heatmapGrob, x = 0, y = 0, height = 5, width = 6)

# Colorbar 
plotGG(plot = heatmapLegendGrob, x = 0.075, y = 5, 
       width = 5.35, height = 0.11)
# Colorbar title
plotText(label = "Relative Expression", fontfamily = "Helvetica",
         fontsize = 8, x = 2.75, y = 5.1, just = "top")

# Age legend
plotRect(x = 5.75, y = 0.15, width = unit(3, "mm"), 
         height = unit(1, "mm"), linecolor = NA, fill = ageColors[6],
         just = "left")
plotRect(x = unit(5.75, "in") + unit(3, "mm"), y = 0.15, width = unit(3, "mm"), 
         height = unit(1, "mm"), linecolor = NA, fill = ageColors[5],
         just = "left")
plotRect(x = unit(5.75, "in") + unit(3*2, "mm"), y = 0.15, width = unit(3, "mm"), 
         height = unit(1, "mm"), linecolor = NA, fill = ageColors[4],
         just = "left")
plotRect(x = unit(5.75, "in") + unit(3*3, "mm"), y = 0.15, width = unit(3, "mm"), 
         height = unit(1, "mm"), linecolor = NA, fill = ageColors[3],
         just = "left")
plotRect(x = unit(5.75, "in") + unit(3*4, "mm"), y = 0.15, width = unit(3, "mm"), 
         height = unit(1, "mm"), linecolor = NA, fill = ageColors[2],
         just = "left")
plotRect(x = unit(5.75, "in") + unit(3*5, "mm"), y = 0.15, width = unit(3, "mm"), 
         height = unit(1, "mm"), linecolor = NA, fill = ageColors[1],
         just = "left")
plotText(label = "30", 
         x = 5.75,
         y = 0.18,
         fontfamily = "Helvetica", fontsize = 4, just = c("left", "top"))
plotText(label = "40", 
         x = unit(5.75, "in") + unit(3, "mm"),
         y = 0.18,
         fontfamily = "Helvetica", fontsize = 4, just = c("left", "top"))
plotText(label = "50", 
         x = unit(5.75, "in") + unit(3*2, "mm"),
         y = 0.18,
         fontfamily = "Helvetica", fontsize = 4, just = c("left", "top"))
plotText(label = "60", 
         x = unit(5.75, "in") + unit(3*3, "mm"),
         y = 0.18,
         fontfamily = "Helvetica", fontsize = 4, just = c("left", "top"))
plotText(label = "70", 
         x = unit(5.75, "in") + unit(3*4, "mm"),
         y = 0.18,
         fontfamily = "Helvetica", fontsize = 4, just = c("left", "top"))
plotText(label = "80", 
         x = unit(5.75, "in") + unit(3*5, "mm"),
         y = 0.18,
         fontfamily = "Helvetica", fontsize = 4, just = c("left", "top"))

# Race Legend
plotRect(x = unit(5.75, "in") + unit(3*6, "mm"), 
         y = 0.35, width = unit(3, "mm"), 
         height = unit(1, "mm"), linecolor = NA, fill = raceColors[5],
         just = "right")
plotRect(x = unit(5.75, "in") + unit(3*5, "mm"), 
         y = 0.35, width = unit(3, "mm"), 
         height = unit(1, "mm"), linecolor = NA, fill = raceColors[4],
         just = "right")
plotRect(x = unit(5.75, "in") + unit(3*4, "mm"), 
         y = 0.35, width = unit(3, "mm"), 
         height = unit(1, "mm"), linecolor = NA, fill = raceColors[3],
         just = "right")
plotRect(x = unit(5.75, "in") + unit(3*3, "mm"), 
         y = 0.35, width = unit(3, "mm"), 
         height = unit(1, "mm"), linecolor = NA, fill = raceColors[2],
         just = "right")
plotRect(x = unit(5.75, "in") + unit(3*2, "mm"), 
         y = 0.35, width = unit(3, "mm"), 
         height = unit(1, "mm"), linecolor = NA, fill = raceColors[1],
         just = "right")
plotText(label = "HI", x = unit(5.75, "in") + unit(3*5.5, "mm"),
         y = 0.425,
         fontsize = 5, fontfamily = "Helvetica")
plotText(label = "C", x = unit(5.75, "in") + unit(3*4.5, "mm"),
         y = 0.425,
         fontsize = 5, fontfamily = "Helvetica")
plotText(label = "BL", x= unit(5.75, "in") + unit(3*3.5, "mm"),
         y = 0.425,
         fontsize = 5, fontfamily = "Helvetica")
plotText(label = "AS", x= unit(5.75, "in") + unit(3*2.5, "mm"),
         y = 0.425,
         fontsize = 5, fontfamily = "Helvetica")
plotText(label = "AR", x= unit(5.75, "in") + unit(3*1.5, "mm"),
         y = 0.425,
         fontsize = 5, fontfamily = "Helvetica")


# Sex legend
plotRect(x = unit(5.75, "in") + unit(3*6, "mm"), 
         y = 0.59, width = unit(3, "mm"), 
         height = unit(1, "mm"), linecolor = NA, fill = sexColors[["M"]],
         just = "right")
plotRect(x = unit(5.75, "in") + unit(3*5, "mm"), 
         y = .59, width = unit(3, "mm"), 
         height = unit(1, "mm"), linecolor = NA, fill = sexColors[["F"]],
         just = "right")
plotText(label = "M", x = unit(5.75, "in") + unit(3*5.5, "mm"),
         y = 0.655,
         fontsize = 5, fontfamily = "Helvetica")
plotText(label = "F", x = unit(5.75, "in") + unit(3*4.5, "mm"),
         y = 0.655,
         fontsize = 5, fontfamily = "Helvetica")


# Condition legend
plotRect(x = unit(5.75, "in") + unit(3*6, "mm"), 
         y = 0.8, width = unit(3, "mm"), 
         height = unit(1, "mm"), linecolor = NA, fill = "#4A4A4A",
         just = "right")
plotRect(x = unit(5.75, "in") + unit(3*5, "mm"), 
         y = 0.8, width = unit(3, "mm"), 
         height = unit(1, "mm"), linecolor = NA, fill = "#B8B8B8",
         just = "right")
plotText(label = "FNF", x = unit(5.75, "in") + unit(3*5.5, "mm"),
         y = 0.85,
         fontsize = 4, fontfamily = "Helvetica")
plotText(label = "CTL", x = unit(5.75, "in") + unit(3*4.5, "mm"),
         y = 0.85,
         fontsize = 4, fontfamily = "Helvetica")

# Add labels of top downregulated and upregulated genes 

plotText(label = upGenes$symbol, x = 5.6, y = seq(0.95, 3.175, length.out = 20),
         fontcolor = "#6B5E27", fontface = "bold",
         fontsize = 6, fontfamily = "Helvetica", just = c("left", "top"))

plotText(label = downGenes$symbol, x = 5.6, y = seq(3.285, 4.8, length.out = 15),
         fontcolor = "#2F4864", fontface = "bold",
         fontsize = 6, fontfamily = "Helvetica", just = c("left", "top"))




dev.off()
