library(DESeq2)
library(ComplexHeatmap)
library(tidyverse)
library(circlize)
library(plotgardener)
library(RColorBrewer)
library(splines)
source("scripts/plotting_utils.R")

# Unique sig untreated and treated genes --------------------------------

ctl_cluster_pval01 <- read_csv("data/ctl_age_pval01clusters.csv") %>%
  pull(gene_id) %>%
  unique()
fnf_cluster_pval01 <- read_csv("data/fnf_age_pval01clusters.csv") %>%
  pull(gene_id) %>%
  unique()


# Check direction of effect in ctl for overlapping and unique fnf  --------
load("data/dds_age_ctl_lrt.rda")

# Overlap
overlappingGenes <- fnf_cluster_pval01[which(fnf_cluster_pval01 %in% ctl_cluster_pval01)]

fnfOverlapping <- read_csv("data/fnf_age_pval01clusters.csv") %>%
  filter(gene_id %in% overlappingGenes) %>%
  distinct(gene_id, .keep_all = TRUE) %>%
  dplyr::select(gene_id, cluster) %>%
  dplyr::rename(fnf_cluster = cluster)

ctlOverlapping <- read_csv("data/ctl_age_pval01clusters.csv") %>%
  filter(gene_id %in% overlappingGenes) %>%
  distinct(gene_id, .keep_all = TRUE) %>%
  dplyr::select(gene_id, cluster) %>%
  dplyr::rename(ctl_cluster = cluster)

ctl_fnf_overlapping <- full_join(ctlOverlapping,
                                 fnfOverlapping,
                                 by = "gene_id")

# Unique
fnf_unique_genes <- fnf_cluster_pval01[which(!fnf_cluster_pval01 %in% ctl_cluster_pval01)]
fnfUnique <- read_csv("data/fnf_age_pval01clusters.csv") %>%
  filter(gene_id %in% fnf_unique_genes) %>%
  distinct(gene_id, .keep_all = TRUE) %>%
  dplyr::select(gene_id, cluster) %>%
  dplyr::rename(fnf_cluster = cluster)

ages <- read_csv("data/ctl_age_pval01clusters.csv") %>%
  pull(Age) %>%
  unique() %>%
  sort()
 

coef_mat <- coef(dds_age_ctl_lrt)
design_mat <- model.matrix(design(dds_age_ctl_lrt), 
                           colData(dds_age_ctl_lrt)) 

for (gene in fnf_unique_genes){
  
  spline_count_data <- plotCounts(dds_age_ctl_lrt, gene, 
                                  intgroup = "Age", returnData = TRUE) %>%
    rownames_to_column(var = "Sample") %>%
    # Build fitted spline from matrices
    mutate(logmu = design_mat %*% coef_mat[gene,],
           gene_id = gene)
  
  # Plot to see pattern
  print(ggplot(spline_count_data, aes(x = Age, y = count)) +
    geom_point(color = "grey50") +
    geom_line(aes(y = logmu), lwd = 0.75, color = "red") +
    scale_y_continuous(name = "Log Expression") +
      ggtitle(gene))
}


# Union of untreated and treated genes ------------------------------------

ctl_cluster_pval01 <- read_csv("data/ctl_age_pval01clusters.csv") %>%
  distinct(gene_id, .keep_all = TRUE) %>%
  dplyr::select(gene_id, cluster)
  

fnf_cluster_pval01 <- read_csv("data/fnf_age_pval01clusters.csv") %>%
  distinct(gene_id, .keep_all = TRUE) %>%
  dplyr::select(gene_id, cluster)

# Join and order by cluster
union_sig_genes <- full_join(ctl_cluster_pval01, fnf_cluster_pval01) %>%
  arrange(desc(cluster))


# Untreated ---------------------------------------------------------------
load("data/dds_age_ctl_lrt.rda")
ages <- read_csv("data/ctl_age_pval01clusters.csv") %>%
  pull(Age) %>%
  unique() %>%
  sort()

# Normalized counts
dds_age_ctl_lrt_norm <- vst(dds_age_ctl_lrt)

# Subset norm counts for significant genes
ctl_cluster_counts <- assay(dds_age_ctl_lrt_norm[union_sig_genes$gene_id,]) %>%
  as.data.frame()
ctl_cluster_counts <- ctl_cluster_counts[match(union_sig_genes$gene_id, 
                                               rownames(ctl_cluster_counts)),]

# Order columns by increasing age
ctl_cluster_counts <- ctl_cluster_counts %>%
  dplyr::select(ends_with(as.character(ages)))

# Scale counts
ctl_cluster_counts_scaled <- t(apply(ctl_cluster_counts, 1, scale))
colnames(ctl_cluster_counts_scaled) <- colnames(ctl_cluster_counts)

# Age, Sex, and Race Clusters
annotations <- as.data.frame(colData(dds_age_ctl_lrt)[,c("Race", "Sex", "Age")])
# Put in same order as matrix
annotations <- annotations[match(colnames(ctl_cluster_counts_scaled), rownames(annotations)),]

age_seqColors <- colorRampPalette(rev(sequential_hcl(n = 9, palette = "Mint")))(length(ages))
names(age_seqColors) <- ages

annotationObjects <- HeatmapAnnotation(
  df = annotations,
  col = list(Race = c("ARAB" = raceColors[1], 
                      "ASIAN" = raceColors[2],
                      "BL" = raceColors[3], 
                      "C" = raceColors[4],
                      "HISP" = raceColors[5],
                      "Unknown" = "grey"),
             Sex = sexColors,
             Age = age_seqColors),
  annotation_name_gp = gpar(fontfamily = "Helvetica",
                            fontsize = 8),
  which = "column"
)


clusters <- HeatmapAnnotation(
  df = union_sig_genes %>%
    column_to_rownames(var = "gene_id"),
  col = list(cluster = c("1" = "#009E8E", 
                         "2" = "#AC7F21")),
  annotation_name_gp = gpar(fontfamily = "Helvetica",
                            fontface = "bold",
                            fontsize = 0),
  which = "row",
  simple_anno_size = unit(3,"mm")
)

age_heatmap_control <- ComplexHeatmap::pheatmap(ctl_cluster_counts_scaled, 
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
                                             fontsize = 8),
                            legend_width = unit(4.55, "in"),
                            grid_height = unit(0.11, "in"),
                            direction = "horizontal"
)

age_heatmap_controlGrob <- grid.grabExpr(draw(age_heatmap_control,
                                              show_annotation_legend = FALSE,
                                              show_heatmap_legend = FALSE,
                                              background = "transparent"))
age_heatmapLegendGrob <- grid.grabExpr(draw(age_heatmapLegend))


# Plot layout together with plotgardener ----------------------------------

pdf(file = "plots/DE_heatmap_age.pdf", width = 5.75, height = 5.25)
pageCreate(width = 5.75, height = 5.25, showGuides = FALSE)


plotGG(plot = age_heatmap_controlGrob, x = 0, y = 0, 
       height = 5, width = 5)


# Colorbar
plotGG(plot = age_heatmapLegendGrob,  x = 0.075, y = 5,
       width = 4.55, height = 0.11)
# Colorbar title
plotText(label = "Relative Expression", fontfamily = "Helvetica",
         fontsize = 8, x = 2.375, y = 5.1, just = "top")


# Race Legend
plotRect(x = unit(5, "in") + unit(3*6, "mm"), 
         y = 0.13, width = unit(3, "mm"), 
         height = unit(1, "mm"), linecolor = NA, fill = raceColors[5],
         just = "right")
plotRect(x = unit(5, "in") + unit(3*5, "mm"), 
         y = 0.13, width = unit(3, "mm"), 
         height = unit(1, "mm"), linecolor = NA, fill = raceColors[4],
         just = "right")
plotRect(x = unit(5, "in") + unit(3*4, "mm"), 
         y = 0.13, width = unit(3, "mm"), 
         height = unit(1, "mm"), linecolor = NA, fill = raceColors[3],
         just = "right")
plotRect(x = unit(5, "in") + unit(3*3, "mm"), 
         y = 0.13, width = unit(3, "mm"), 
         height = unit(1, "mm"), linecolor = NA, fill = raceColors[2],
         just = "right")
plotRect(x = unit(5, "in") + unit(3*2, "mm"), 
         y = 0.13, width = unit(3, "mm"), 
         height = unit(1, "mm"), linecolor = NA, fill = raceColors[1],
         just = "right")
plotText(label = "HI", x = unit(5, "in") + unit(3*5.5, "mm"),
         y = 0.2,
         fontsize = 5, fontfamily = "Helvetica")
plotText(label = "C", x = unit(5, "in") + unit(3*4.5, "mm"),
         y = 0.2,
         fontsize = 5, fontfamily = "Helvetica")
plotText(label = "BL", x= unit(5, "in") + unit(3*3.5, "mm"),
         y = 0.2,
         fontsize = 5, fontfamily = "Helvetica")
plotText(label = "AS", x= unit(5, "in") + unit(3*2.5, "mm"),
         y = 0.2,
         fontsize = 5, fontfamily = "Helvetica")
plotText(label = "AR", x= unit(5, "in") + unit(3*1.5, "mm"),
         y = 0.2,
         fontsize = 5, fontfamily = "Helvetica")

# Sex legend
plotRect(x = unit(5, "in") + unit(3*6, "mm"), 
         y = 0.35, width = unit(3, "mm"), 
         height = unit(1, "mm"), linecolor = NA, fill = sexColors[["M"]],
         just = "right")
plotRect(x = unit(5, "in") + unit(3*5, "mm"), 
         y = .35, width = unit(3, "mm"), 
         height = unit(1, "mm"), linecolor = NA, fill = sexColors[["F"]],
         just = "right")
plotText(label = "M", x = unit(5, "in") + unit(3*5.5, "mm"),
         y = 0.41,
         fontsize = 5, fontfamily = "Helvetica")
plotText(label = "F", x = unit(5, "in") + unit(3*4.5, "mm"),
         y = 0.41,
         fontsize = 5, fontfamily = "Helvetica")

# Age legend
xcoord <- unit(5, "in")
for (i in 1:length(age_seqColors)){
  plotRect(x = xcoord, y = 0.59, width = unit(18/38, "mm"),
           height = unit(1, "mm"), linecolor = NA,
           fill = age_seqColors[i], just = "left")
  xcoord <- xcoord + unit(18/38, "mm")
}

plotText(label = "34", 
         x = 5,
         y = 0.625,
         fontfamily = "Helvetica", fontsize = 5, just = c("left", "top"))

plotText(label = "84", 
         x = unit(5, "in") + unit(3*6, "mm"),
         y = 0.625,
         fontfamily = "Helvetica", fontsize = 5, just = c("right", "top"))
dev.off()

