library(DESeq2)
library(ComplexHeatmap)
library(tidyverse)
library(circlize)
library(plotgardener)
library(RColorBrewer)
library(splines)
library(rrvgo)
library(ggvenn)
source("scripts/plotting_utils.R")


# Functions ---------------------------------------------------------------

splineFitAge <- function(gene, coef_mat, design_mat){
  
  gene_spline_fit <- data.frame("Sample" = rownames(design_mat),
                                "logmu" = design_mat %*% coef_mat[gene,],
                                "gene_id" = gene) %>%
    separate_wider_delim(cols = "Sample", delim = "_", 
                         names = c(NA, NA, NA, NA, NA, NA, "Age"),
                         cols_remove = FALSE) %>%
    mutate(Age = as.numeric(Age))
  
  return(gene_spline_fit)
  
}

gene_enrichment_permTest <- function(group, nPerm = 1000){
  
  gene_dataset <- switch(group,
                         fnf = fnf_de_genes,
                         oa = oa_de_genes,
                         fnf_oa = fnf_oa_de_genes)
  
  age_dataset <- switch(group,
                        fnf = age_fnf_degenes,
                        oa = age_oa_degenes,
                        fnf_oa = age_fnfoa_degenes)
  
  randomOverlaps <- c() 
  for (i in 1:nPerm){
    # Randomly select number of sex_degenes_union from all genes
    random_genes <- de_genes_results[sample(nrow(de_genes_results),
                                            nrow(union_sig_genes)),]
    
    # Check how many overlap with group
    groupOverlap <- length(which(random_genes$symbol %in% gene_dataset$symbol))
    
    # Append
    randomOverlaps[i] <- groupOverlap
  }
  
  # Calculate p-val
  pval <- length(which(randomOverlaps > nrow(age_dataset)))/nPerm
  return(list("pval" = pval, "randomOverlaps" = randomOverlaps))
  }
  

get_gene_age_Counts <- function(gene, dds){
  
  # Get coef matrix and design matrix for fitted spline
  coef_mat <- coef(dds)
  design_mat <- model.matrix(design(dds), colData(dds))
  
  geneCounts <- plotCounts(dds, gene = gene, 
                           intgroup = "Age", 
                           normalized = TRUE,
                           returnData = TRUE) %>%
    remove_rownames() %>%
    mutate(gene_id = gene, 
           logmu = design_mat %*% coef_mat[gene,])
  
  return(geneCounts)
}

# HEATMAP -----------------------------------------------------------------

# Union of untreated and treated genes
ctl_cluster_pval01 <- read_csv("data/ctl_age_pval01clusters.csv") %>%
  distinct(gene_id, .keep_all = TRUE) %>%
  dplyr::select(gene_id, cluster)

fnf_cluster_pval01 <- read_csv("data/fnf_age_pval01clusters.csv") %>%
  distinct(gene_id, .keep_all = TRUE) %>%
  dplyr::select(gene_id, cluster)

# Join and order by cluster
union_sig_genes <- full_join(ctl_cluster_pval01, fnf_cluster_pval01) %>%
  arrange(desc(cluster))

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
  col = list(cluster = ageClusterColors),
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


# Plot layout together with plotgardener
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
plotRect(x = unit(5, "in") + unit(3*4, "mm"), 
         y = 0.13, width = unit(3, "mm"), 
         height = unit(1, "mm"), linecolor = NA, fill = raceColors[5],
         just = "left")
plotRect(x = unit(5, "in") + unit(3*3, "mm"), 
         y = 0.13, width = unit(3, "mm"), 
         height = unit(1, "mm"), linecolor = NA, fill = raceColors[4],
         just = "left")
plotRect(x = unit(5, "in") + unit(3*2, "mm"), 
         y = 0.13, width = unit(3, "mm"), 
         height = unit(1, "mm"), linecolor = NA, fill = raceColors[3],
         just = "left")
plotRect(x = unit(5, "in") + unit(3*1, "mm"), 
         y = 0.13, width = unit(3, "mm"), 
         height = unit(1, "mm"), linecolor = NA, fill = raceColors[2],
         just = "left")
plotRect(x = unit(5, "in"), 
         y = 0.13, width = unit(3, "mm"), 
         height = unit(1, "mm"), linecolor = NA, fill = raceColors[1],
         just = "left")
plotText(label = "HI", x = unit(5, "in") + unit(3*4.5, "mm"),
         y = 0.2,
         fontsize = 5, fontfamily = "Helvetica")
plotText(label = "C", x = unit(5, "in") + unit(3*3.5, "mm"),
         y = 0.2,
         fontsize = 5, fontfamily = "Helvetica")
plotText(label = "BL", x= unit(5, "in") + unit(3*2.5, "mm"),
         y = 0.2,
         fontsize = 5, fontfamily = "Helvetica")
plotText(label = "AS", x= unit(5, "in") + unit(3*1.5, "mm"),
         y = 0.2,
         fontsize = 5, fontfamily = "Helvetica")
plotText(label = "AR", x= unit(5, "in") + unit(3*0.5, "mm"),
         y = 0.2,
         fontsize = 5, fontfamily = "Helvetica")

# Sex legend
plotRect(x = unit(5, "in") + unit(3, "mm"), 
         y = 0.35, width = unit(3, "mm"), 
         height = unit(1, "mm"), linecolor = NA, fill = sexColors[["M"]],
         just = "left")
plotRect(x = unit(5, "in"), 
         y = .35, width = unit(3, "mm"), 
         height = unit(1, "mm"), linecolor = NA, fill = sexColors[["F"]],
         just = "left")
plotText(label = "M", x = unit(5, "in") + unit(3*1.5, "mm"),
         y = 0.41,
         fontsize = 5, fontfamily = "Helvetica")
plotText(label = "F", x = unit(5, "in") + unit(3*0.5, "mm"),
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

# CLUSTER LINE PLOTS ------------------------------------------------------

ctl_cluster_pval01 <- read_csv("data/ctl_age_pval01clusters.csv") %>%
  distinct(gene_id, .keep_all = TRUE) %>%
  dplyr::select(gene_id, cluster)

fnf_cluster_pval01 <- read_csv("data/fnf_age_pval01clusters.csv") %>%
  distinct(gene_id, .keep_all = TRUE) %>%
  dplyr::select(gene_id, cluster)

# Join and order by cluster
union_sig_genes <- full_join(ctl_cluster_pval01, fnf_cluster_pval01) %>%
  arrange(desc(cluster))

# Get spline fit data from control 
load("data/dds_age_ctl_lrt.rda")

coef_mat <- coef(dds_age_ctl_lrt)
design_mat <- model.matrix(design(dds_age_ctl_lrt), colData(dds_age_ctl_lrt))

ctl_sig_gene_fits <- bind_rows(lapply(union_sig_genes$gene_id, 
                                      splineFitAge, coef_mat, design_mat)) %>%
  left_join(union_sig_genes, by = "gene_id")

# Full cluster line plots
for (clust in 1:2){
  
  if (clust == 1){
    median_color <- ageClusterColors[["1"]]
    limits <- c(0, 16)
    breaks <- seq(4, 16, 4)
    
  } else if (clust == 2){
    median_color <- ageClusterColors[["2"]]
    limits <- c(0,14)
    breaks <- seq(2, 14, 4) 
  } 
  
  ctl_cluster <- ctl_sig_gene_fits %>%
    filter(cluster == clust)
  
  ggplot(ctl_cluster, aes(x = Age, y = logmu, color = gene_id)) +
    geom_line(lwd = 0.5) +
    scale_color_manual(values = rep("grey75", nrow(ctl_cluster))) +
    stat_summary(geom = "line", fun = "median", color = median_color, lwd = 1) +
    geom_vline(xintercept = 62, lty = 2) +
    scale_x_continuous(breaks = c(40, 50, 60, 70, 80)) +
    scale_y_continuous(name = "Log Expression Fit",
                       limits = limits,
                       breaks = breaks) +
    coord_cartesian(clip = "off") +
    theme(panel.background = element_blank(),
          text = element_text(family = "Helvetica"),
          panel.border = element_rect(fill = NA),
          legend.position = "none",
          axis.text = element_text(color = "black", size = 12),
          axis.title = element_text(size = 12),
          axis.line = element_line(color = "black"),
          axis.ticks.x = element_blank(),
          axis.ticks.length.y = unit(-0.1, "cm"))
  
  ggsave(filename = paste0("plots/splineFit_unionsigGenes_cluster", clust, "_splines.pdf"),
         width = 6, height = 6, units = "in")
  
}

# GO AND KEGG BARPLOTS FOR EACH CLUSTER -----------------------------------

## GO terms
age_cluster1_go <- 
  read_delim("data/homer_age_sig_cluster1/biological_process.txt") %>%
  mutate(pval = exp(1)^logP) %>%
  filter(pval < 0.01)
age_cluster1_go <- reduceGO(age_cluster1_go,
                     category = "clust1")

age_cluster2_go <-
  read_delim("data/homer_age_sig_cluster2/biological_process.txt") %>%
  mutate(pval = exp(1)^logP) %>%
  filter(pval < 0.01)
age_cluster2_go <- reduceGO(age_cluster2_go,
                            category = "clust2")

age_go_plotting <- bind_rows(age_cluster1_go, age_cluster2_go)
age_go_plotting$parentTerm <- factor(age_go_plotting$parentTerm, levels = rev(age_go_plotting$parentTerm))
age_go_plotting$category <- factor(age_go_plotting$category, levels = c("clust2", "clust1"))

# Plot all in barplot
ggplot(age_go_plotting, aes(x = `-log10pval`, y = parentTerm, fill = category)) +
  geom_vline(xintercept = 5, color = "grey75", alpha = 0.4) +
  geom_vline(xintercept = 10, color = "grey75", alpha = 0.4) +
  geom_vline(xintercept = 15, color = "grey75", alpha = 0.4) +
  geom_bar(stat = "identity") +
  scale_x_continuous(limits = c(0, 15), expand = c(0, 0), name = "-log~10~pval",
                     breaks = seq(0, 15, 5)) +
  scale_fill_manual(values = c(ageClusterColors[["2"]], ageClusterColors[["1"]])) +
  coord_cartesian(clip = "off") +
  facet_wrap(~category, ncol = 1, strip.position = "left", scales = "free_y") +
  geom_text(aes(x = 0, label = parentTerm), hjust = 0, family = "Helvetica") +
  theme(panel.background = element_blank(),
        text = element_text(family = "Helvetica"),
        legend.position = "None",
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_markdown(),
        axis.text.x = element_text(color = "black", size = 10),
        strip.background = element_blank(),
        strip.text.y.left = element_blank(),
        axis.ticks = element_blank(),
        axis.line.x = element_line(),
        strip.text = element_text(size = 14, color = "black"),
        panel.spacing = unit(0, "mm"), 
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14)) +
  ggtitle("GO Terms")

ggsave(filename = "plots/age_GO_barplots.pdf", width = 6, height = 10, units = "in")

### KEGG PATHWAYS
clust1_kegg_data <- read_delim("data/homer_age_sig_cluster1/kegg.txt") %>%
  mutate(pval = exp(1)^logP) %>%
  filter(pval < 0.01) %>%
  distinct(Term, .keep_all = TRUE) %>%
  mutate(`-log10pval` = -log10(pval)) %>%
  mutate(category = "clust1")
clust2_kegg_data <- read_delim("data/homer_age_sig_cluster2/kegg.txt") %>%
  mutate(pval = exp(1)^logP) %>%
  filter(pval < 0.01) %>%
  distinct(Term, .keep_all = TRUE) %>%
  mutate(`-log10pval` = -log10(pval)) %>%
  mutate(category = "clust2")

age_kegg_plotting <- bind_rows(clust1_kegg_data, clust2_kegg_data)
age_kegg_plotting$Term <- factor(age_kegg_plotting$Term, levels = age_kegg_plotting$Term)
age_kegg_plotting$category <- factor(age_kegg_plotting$category, levels = c("clust2", "clust1"))

ggplot(age_kegg_plotting, aes(x = `-log10pval`, y = Term, fill = category)) +
  geom_vline(xintercept = 1, color = "grey75", alpha = 0.4) +
  geom_vline(xintercept = 2, color = "grey75", alpha = 0.4) +
  geom_vline(xintercept = 3, color = "grey75", alpha = 0.4) +
  geom_vline(xintercept = 4, color = "grey75", alpha = 0.4) +
  geom_vline(xintercept = 5, color = "grey75", alpha = 0.4) +
  geom_bar(stat = "identity") +
  scale_x_continuous(expand = c(0, 0), name = "-log~10~pval", limits = c(0, 5),
                     breaks = seq(0, 5, 1)) +
  scale_fill_manual(values = c(ageClusterColors[["2"]], ageClusterColors[["1"]])) +
  facet_wrap(~category, ncol = 1, strip.position = "left", scales = "free_y") +
  geom_text(aes(x = 0, label = Term), hjust = 0, family = "Helvetica") +
  theme(panel.background = element_blank(),
        text = element_text(family = "Helvetica"),
        legend.position = "None",
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_markdown(),
        axis.text.x = element_text(color = "black", size = 10),
        strip.background = element_blank(),
        axis.ticks = element_blank(),
        axis.line.x = element_line(),
        strip.text = element_blank(),
        panel.spacing = unit(0, "mm"), 
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14)) +
  ggtitle("KEGG Pathways")

ggsave(filename = "plots/age_KEGG_barplots.pdf", width = 4, height = 6, units = "in")

# VENN DIAGRAM OF PBS AND FNF AGE GENES -----------------------------------
ctl_cluster_pval01 <- read_csv("data/ctl_age_pval01clusters.csv") %>%
  distinct(gene_id, .keep_all = TRUE) %>%
  dplyr::select(gene_id, cluster)

fnf_cluster_pval01 <- read_csv("data/fnf_age_pval01clusters.csv") %>%
  distinct(gene_id, .keep_all = TRUE) %>%
  dplyr::select(gene_id, cluster)


pbs_fnf_age_genes <- tibble(values = unique(c(ctl_cluster_pval01$gene_id, fnf_cluster_pval01$gene_id))) %>%
  mutate(PBS = values %in% ctl_cluster_pval01$gene_id,
         FNF = values %in% fnf_cluster_pval01$gene_id)


venn_diagram <- ggplot(pbs_fnf_age_genes, aes(A = PBS, B = FNF)) +
  geom_venn(set_names = c("PBS", "FN-f"), 
            fill_color = c(log2fcColors[["-"]], log2fcColors[["+"]]), 
            stroke_color = NA, auto_scale = TRUE, show_percentage = FALSE,
            text_size = 5, set_name_size = 6) +
  coord_fixed()  +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank())
venn_font(venn_diagram, font = "Helvetica")

ggsave(filename = "plots/age_pbs_fnf_venndiagram.pdf", width = 4, height = 4, units = "in")



# OVERLAP WITH FNF AND OA -------------------------------------------------

ctl_cluster_pval01 <- read_csv("data/ctl_age_pval01clusters.csv") %>%
  distinct(gene_id, .keep_all = TRUE) %>%
  dplyr::select(gene_id, symbol, cluster) %>%
  dplyr::rename(ctl_cluster = cluster)

fnf_cluster_pval01 <- read_csv("data/fnf_age_pval01clusters.csv") %>%
  distinct(gene_id, .keep_all = TRUE) %>%
  dplyr::select(gene_id, symbol, cluster) %>%
  dplyr::rename(fnf_cluster = cluster)

union_sig_genes <- full_join(ctl_cluster_pval01, fnf_cluster_pval01, 
                             by = c("gene_id", "symbol"))

fnf_de_genes <- read_csv("data/sig_deGenes_pval01_l2fc1.csv")
oa_de_genes <- read_csv("data/RAAK_genes.csv",
                        col_select = c("ENSEMBL", "HGNC", "RAAK_PVAL",
                                       "RAAK_FC", "RAAK_LFC")) %>%
  dplyr::rename(symbol = HGNC) %>%
  dplyr::rename(log2FoldChange = RAAK_LFC) %>%
  dplyr::rename(gene_id = ENSEMBL)
fnf_oa_de_genes <- fnf_de_genes %>%
  filter(symbol %in% oa_de_genes$symbol)

# Get age de genes only differential in FN-f
age_fnf_degenes <- union_sig_genes %>%
  filter(symbol %in% fnf_de_genes$symbol & !gene_id %in% oa_de_genes$symbol)

# All age de genes differential in FN-f
age_fnf_degenes_all <- union_sig_genes %>%
  filter(symbol %in% fnf_de_genes$symbol) 

# Get age de genes only differential in OA
age_oa_degenes <- union_sig_genes %>%
  filter(!symbol %in% fnf_de_genes$symbol & symbol %in% oa_de_genes$symbol) 

# All age de genes differential in OA
age_oa_degenes_all <- union_sig_genes %>%
  filter(symbol %in% oa_de_genes$symbol) %>%
  left_join(oa_de_genes %>% dplyr::select(symbol, log2FoldChange)) %>%
  dplyr::rename(oa_log2FoldChange = log2FoldChange) %>%
  mutate(oa_group = ifelse(oa_log2FoldChange < 0, "down", "up"))

# Get age de genes differential in FN-f and OA
age_fnfoa_degenes <- union_sig_genes %>%
  filter(symbol %in% fnf_de_genes$symbol & symbol %in% oa_de_genes$symbol)

## Permutation tests for enrichment
de_genes_results <- read_csv("data/de_genes_results.csv")

fnf_enrichment <- gene_enrichment_permTest(group = "fnf")
oa_enrichment <- gene_enrichment_permTest(group = "oa")
fnf_oa_enrichment <- gene_enrichment_permTest(group = "fnf_oa")

overlapData <- data.frame(group = c("FN-f", "FN-f and OA", "OA"),
                          number = c(nrow(age_fnf_degenes),
                                     nrow(age_fnfoa_degenes),
                                     nrow(age_oa_degenes)),
                          pval = c(fnf_enrichment$pval, fnf_oa_enrichment$pval, oa_enrichment$pval),
                          expected_value = c(median(fnf_enrichment$randomOverlaps),
                                             median(fnf_oa_enrichment$randomOverlaps),
                                             median(oa_enrichment$randomOverlaps)),
                          expected_value_q1 = c(quantile(fnf_enrichment$randomOverlaps, probs = 0.25),
                                                quantile(fnf_oa_enrichment$randomOverlaps, probs = 0.25),
                                                quantile(oa_enrichment$randomOverlaps, probs = 0.25)),
                          expected_value_q3 = c(quantile(fnf_enrichment$randomOverlaps, probs = 0.75),
                                                quantile(fnf_oa_enrichment$randomOverlaps, probs = 0.75),
                                                quantile(oa_enrichment$randomOverlaps, probs = 0.75)))

### Line plots of overlapping OA genes
load("data/dds_age_ctl_lrt.rda")
load("data/dds_age_fnf_lrt.rda")

all_gene_age_counts <- list()
for (geneRow in 1:nrow(age_oa_degenes_all)){
  gene <- age_oa_degenes_all[geneRow, ]
  
  ctl_gene_counts <- get_gene_age_Counts(gene = gene[["gene_id"]], dds = dds_age_ctl_lrt) %>%
    mutate(condition = "PBS")
  fnf_gene_counts <- get_gene_age_Counts(gene = gene[["gene_id"]], dds = dds_age_fnf_lrt) %>%
    mutate(condition = "FN-f")
  
  all_gene_counts <- bind_rows(ctl_gene_counts, fnf_gene_counts)
  all_gene_age_counts[[gene[["gene_id"]]]] <- all_gene_counts
}

all_gene_age_counts <- left_join(bind_rows(all_gene_age_counts),
                                 age_oa_degenes_all %>%
                                   dplyr::select(symbol, gene_id, oa_group),
                                 by = "gene_id")
all_gene_age_counts$symbol <- factor(all_gene_age_counts$symbol,
                                     levels = c("CXCL14", "PAPPA", "GPC5"))

all_gene_age_counts$oa_group <- factor(all_gene_age_counts$oa_group,
                                       levels = c("up", "down"))
all_gene_age_counts$condition <- factor(all_gene_age_counts$condition,
                                        levels = c("PBS", "FN-f"))


sig_age_genes <- age_oa_degenes_all %>%
  pivot_longer(cols = ends_with("cluster"), 
               names_to = "condition", 
               values_to = "cluster") %>%
  filter(!is.na(cluster)) %>%
  mutate(condition = gsub("_cluster", "", condition)) %>%
  mutate(condition = ifelse(condition == "ctl", "PBS", "FN-f")) %>%
  mutate(cluster = as.character(cluster))
sig_age_genes$condition <- factor(sig_age_genes$condition,
                                  levels = c("PBS", "FN-f"))



up_genes <- ggplot(all_gene_age_counts %>% 
         filter(oa_group == "up"), aes(x = Age, y = log2(count))) +
  geom_point(color = "grey25") +
  #stat_poly_line(se = FALSE) +
  #stat_poly_eq() +
  geom_line(aes(x = Age, y = logmu), color = ageClusterColors[["2"]], lwd = 1) +
  geom_star(data = sig_age_genes %>% 
              filter(oa_group == "up"),
            aes(x = 60, y = Inf, starshape = cluster),
            fill = ageClusterColors[["2"]],
            size = 3) +
  ggh4x::facet_grid2(condition ~ symbol, switch = "both", scales = "free",
                     independent = "all") +
  coord_cartesian(clip = "off") +
  scale_y_continuous(breaks = scales::breaks_pretty(3), 
                     name = "log2(normalized counts)") +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        text = element_text(family = "Helvetica"),
        strip.placement = "outside",
        axis.line = element_line(linewidth = 0.25),
        axis.ticks.length = unit(-0.1, "cm"),
        axis.title = element_text(size = 8),
        axis.text = element_text(color = "black"),
        strip.background = element_blank(),
        strip.text.y.left = element_text(angle = 360, size = 10),
        strip.text.x.bottom = element_text(face = "bold", size = 12),
        axis.ticks = element_line(color = "black", linewidth = 0.25),
        panel.spacing.y = unit(0.5, "cm"),
        plot.title = element_text(hjust = 0.5),
        legend.position = "None") +
  ggtitle("Up in OA")

down_genes <- ggplot(all_gene_age_counts %>% 
                     filter(oa_group == "down"), aes(x = Age, y = log2(count))) +
  geom_point(color = "grey25") +
  #stat_poly_line(se = FALSE) +
  #stat_poly_eq() +
  geom_line(aes(x = Age, y = logmu), color = ageClusterColors[["2"]], lwd = 1) +
  geom_star(data = sig_age_genes %>% 
              filter(oa_group == "down"),
            aes(x = 60, y = Inf, starshape = cluster),
            fill = ageClusterColors[["2"]],
            size = 3) +
  ggh4x::facet_grid2(condition ~ symbol, switch = "both", scales = "free",
                     independent = "all") +
  coord_cartesian(clip = "off") +
  scale_y_continuous(breaks = scales::breaks_pretty(3), 
                     name = "log2(normalized counts)") +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        text = element_text(family = "Helvetica"),
        axis.title.y = element_blank(),
        strip.placement = "outside",
        axis.line = element_line(linewidth = 0.25),
        axis.ticks.length = unit(-0.1, "cm"),
        axis.text = element_text(color = "black"),
        strip.background = element_blank(),
        strip.text.y.left = element_blank(),
        strip.text.x.bottom = element_text(face = "bold", size = 12),
        axis.ticks = element_line(color = "black", linewidth = 0.25),
        panel.spacing.y = unit(0.5, "cm"),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        legend.position = "None") +
  ggtitle("Down in OA")

up_genes + down_genes + plot_layout(widths = c(2, 1))

ggsave(filename = "plots/oa_age_genes_counts.pdf", 
       width = 8, height = 5, units = "in")

## Splitting into age groups and making box plots

gene_maxCounts <- all_gene_age_counts %>%
  group_by(symbol) %>%
  summarise(maxCount = max(log2(count)))


all_gene_age_counts <- all_gene_age_counts %>%
  mutate(Age_group = ifelse(Age < 50, "<50", ifelse(Age < 65, "50-65", "65+")))
  
all_gene_age_counts$Age_group <- factor(all_gene_age_counts$Age_group,
                                        levels = c("<50", "50-65", "65+"))

sig_age_genes <- left_join(sig_age_genes, gene_maxCounts, by = "symbol")
sig_age_genes <- merge(sig_age_genes, data.frame("Age_group" = c("<50", "50-65", "65+")))


up_age_boxplots <- ggplot(all_gene_age_counts %>%
         filter(oa_group == "up"), aes(x = Age_group, y = log2(count))) +
  geom_boxplot(outlier.shape = NA, fill = ageClusterColors[["2"]]) +
  geom_line(data = sig_age_genes %>% 
              filter(oa_group == "up"), 
            aes(x = Age_group, y = maxCount + 0.5, group = symbol), inherit.aes = FALSE,
            linewidth = 0.5) +
  geom_star(data = sig_age_genes %>% 
              filter(oa_group == "up"),
            aes(x = 2, y = Inf, starshape = cluster),
            fill = ageClusterColors[["2"]],
            size = 3) +
  ggh4x::facet_grid2(condition ~ symbol, switch = "both", scales = "free",
                     independent = "all") +
  scale_y_continuous(breaks = scales::breaks_pretty(3), 
                     name = "log2(normalized counts)") +
  coord_cartesian(clip = "off") +
  theme(strip.placement = "outside",
        axis.line = element_line(linewidth = 0.25),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(color = "black", linewidth = 0.25),
        axis.ticks.length.y = unit(-0.1, "cm"),
        axis.title.y = element_text(vjust = -14, size = 8, family= "Helvetica"),
        text = element_text(family = "Helvetica"),
        axis.text = element_text(color = "black"),
        strip.background = element_blank(),
        strip.text.y.left = element_text(angle = 360, size = 10),
        strip.text.x.bottom = element_text(face = "bold", size = 12),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none",
        panel.spacing.y = unit(0.5, "cm"),
        plot.title = element_text(hjust = 0.5)) +
  ggtitle("Up in OA")

down_age_boxplots <- ggplot(all_gene_age_counts %>%
         filter(oa_group == "down"), aes(x = Age_group, y = log2(count))) +
  geom_boxplot(outlier.shape = NA, fill = ageClusterColors[["2"]]) +
  geom_line(data = sig_age_genes %>% 
              filter(oa_group == "down"), 
            aes(x = Age_group, y = maxCount + 0.5, group = symbol), inherit.aes = FALSE,
            linewidth = 0.5) +
  geom_star(data = sig_age_genes %>% 
              filter(oa_group == "down"),
            aes(x = 2, y = Inf, starshape = cluster),
            fill = ageClusterColors[["2"]],
            size = 3) +
  ggh4x::facet_grid2(condition ~ symbol, switch = "both", scales = "free",
                     independent = "all") +
  scale_y_continuous(breaks = scales::breaks_pretty(3), 
                     name = "log2(normalized counts)") +
  coord_cartesian(clip = "off") +
  theme(strip.placement = "outside",
        axis.line = element_line(linewidth = 0.25),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(color = "black", linewidth = 0.25),
        axis.ticks.length.y = unit(-0.1, "cm"),
        axis.title.y = element_blank(),
        text = element_text(family = "Helvetica"),
        axis.text = element_text(color = "black"),
        strip.background = element_blank(),
        strip.text.y.left =element_blank(),
        strip.text.x.bottom = element_text(face = "bold", size = 12),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none",
        panel.spacing.y = unit(0.5, "cm"),
        plot.title = element_text(hjust = 0.5)) +
  ggtitle("Down in OA")

up_age_boxplots + down_age_boxplots + plot_layout(widths = c(2, 1))
ggsave(filename = "plots/oa_age_genes_boxplots.pdf", 
       width = 9, height = 6, units = "in")
