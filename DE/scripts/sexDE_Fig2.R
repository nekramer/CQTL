library(DESeq2)
library(ComplexHeatmap)
library(tidyverse)
library(circlize)
library(patchwork)
library(ggtext)
library(ggstar)
library(plotgardener)
library(ggvenn)
source("scripts/plotting_utils.R")

# Functions ---------------------------------------------------------------

get_gene_sex_Counts <- function(gene, dds){
  
  geneCounts <- plotCounts(dds, gene = gene, 
                           intgroup = "Sex", 
                           normalized = TRUE,
                           returnData = TRUE) %>%
    remove_rownames() %>%
    mutate(gene_id = gene)
  
  return(geneCounts)
}

gene_enrichment_permTest <- function(group, nPerm = 1000){
  
  gene_dataset <- switch(group,
                    fnf = fnf_de_genes,
                    oa = oa_de_genes,
                    fnf_oa = fnf_oa_de_genes)
  
  sex_dataset <- switch(group,
                         fnf = sex_fnf_degenes,
                         oa = sex_oa_degenes,
                         fnf_oa = sex_fnfoa_degenes)
  
  randomOverlaps <- c() 
  for (i in 1:nPerm){
    # Randomly select number of sex_degenes_union from all genes
    random_genes <- de_genes_results[sample(nrow(de_genes_results),
                                            nrow(sex_degenes_union)),]
    
    # Check how many overlap with group
    groupOverlap <- length(which(random_genes$symbol %in% gene_dataset$symbol))
    
    # Append
    randomOverlaps[i] <- groupOverlap
  }
 
  # Calculate p-val
  pval <- length(which(randomOverlaps >= nrow(sex_dataset)))/nPerm
  return(pval)
  
}

# HEATMAP -----------------------------------------------------------------

ctl_sig_genes <- read_csv("data/sex_de/ctl_sexDE_pval01.csv")
fnf_sig_genes <- read_csv("data/sex_de/fnf_sexDE_pval01.csv")

union_sig_genes <- union(ctl_sig_genes$gene_id, fnf_sig_genes$gene_id)

load("data/sex_de/dds_sex_ctl.rda")

# Normalized counts
dds_sex_ctl_norm <- vst(dds_sex_ctl)

## Read in data from ctl selecting union sig genes and order by log2FC
untreated_sex_degenes <- read_csv("data/sex_de/ctl_sex_shrink.csv") %>%
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

annotationObjects <- ComplexHeatmap::HeatmapAnnotation(
  df = annotations,
  col = list(Race = c("ARAB" = raceColors[1], 
                      "ASIAN" = raceColors[2],
                      "BL" = raceColors[3], 
                      "C" = raceColors[4],
                      "HISP" = raceColors[5],
                      "Unknown" = "grey"),
             Age_group = c("31-40" = ageColors[1],
                           "41-50" = ageColors[2],
                           "51-60" = ageColors[3],
                           "61-70" = ageColors[4],
                           "71-80" = ageColors[5],
                           "81-90" = ageColors[6]),
             Sex = sexColors),
  annotation_name_gp = gpar(fontfamily = "Helvetica",
                            fontsize = 8),
  annotation_label = c("Age", "Race", "Sex"),
  which = "column"
)

cluster_annotations <- untreated_sex_degenes %>% 
  dplyr::select(gene_id, log2FC_dir) %>%
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


pdf(file = "plots/sexDE_Fig2/DE_heatmap_sex.pdf", width = 5.75, height = 5.25)
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
plotText(label = "31", 
         x = 5,
         y = 0.18,
         fontfamily = "Helvetica", fontsize = 4, just = c("left", "top"))
plotText(label = "41", 
         x = unit(5, "in") + unit(3, "mm"),
         y = 0.18,
         fontfamily = "Helvetica", fontsize = 4, just = c("left", "top"))
plotText(label = "51", 
         x = unit(5, "in") + unit(3*2, "mm"),
         y = 0.18,
         fontfamily = "Helvetica", fontsize = 4, just = c("left", "top"))
plotText(label = "61", 
         x = unit(5, "in") + unit(3*3, "mm"),
         y = 0.18,
         fontfamily = "Helvetica", fontsize = 4, just = c("left", "top"))
plotText(label = "71", 
         x = unit(5, "in") + unit(3*4, "mm"),
         y = 0.18,
         fontfamily = "Helvetica", fontsize = 4, just = c("left", "top"))
plotText(label = "81", 
         x = unit(5, "in") + unit(3*5, "mm"),
         y = 0.18,
         fontfamily = "Helvetica", fontsize = 4, just = c("left", "top"))

# Race Legend
plotRect(x = unit(5, "in") + unit(3*4, "mm"), 
         y = 0.35, width = unit(3, "mm"), 
         height = unit(1, "mm"), linecolor = NA, fill = raceColors[5],
         just = "left")
plotRect(x = unit(5, "in") + unit(3*3, "mm"), 
         y = 0.35, width = unit(3, "mm"), 
         height = unit(1, "mm"), linecolor = NA, fill = raceColors[4],
         just = "left")
plotRect(x = unit(5, "in") + unit(3*2, "mm"), 
         y = 0.35, width = unit(3, "mm"), 
         height = unit(1, "mm"), linecolor = NA, fill = raceColors[3],
         just = "left")
plotRect(x = unit(5, "in") + unit(3*1, "mm"), 
         y = 0.35, width = unit(3, "mm"), 
         height = unit(1, "mm"), linecolor = NA, fill = raceColors[2],
         just = "left")
plotRect(x = unit(5, "in"), 
         y = 0.35, width = unit(3, "mm"), 
         height = unit(1, "mm"), linecolor = NA, fill = raceColors[1],
         just = "left")
plotText(label = "HI", x = unit(5, "in") + unit(3*4.5, "mm"),
         y = 0.425,
         fontsize = 5, fontfamily = "Helvetica")
plotText(label = "C", x = unit(5, "in") + unit(3*3.5, "mm"),
         y = 0.425,
         fontsize = 5, fontfamily = "Helvetica")
plotText(label = "BL", x= unit(5, "in") + unit(3*2.5, "mm"),
         y = 0.425,
         fontsize = 5, fontfamily = "Helvetica")
plotText(label = "AS", x= unit(5, "in") + unit(3*1.5, "mm"),
         y = 0.425,
         fontsize = 5, fontfamily = "Helvetica")
plotText(label = "AR", x= unit(5, "in") + unit(3*0.5, "mm"),
         y = 0.425,
         fontsize = 5, fontfamily = "Helvetica")


# Sex legend
plotRect(x = unit(5, "in") + unit(3, "mm"), 
         y = 0.59, width = unit(3, "mm"), 
         height = unit(1, "mm"), linecolor = NA, fill = sexColors[["M"]],
         just = "left")
plotRect(x = unit(5, "in"), 
         y = .59, width = unit(3, "mm"), 
         height = unit(1, "mm"), linecolor = NA, fill = sexColors[["F"]],
         just = "left")
plotText(label = "M", x = unit(5, "in") + unit(3*1.5, "mm"),
         y = 0.655,
         fontsize = 5, fontfamily = "Helvetica")
plotText(label = "F", x = unit(5, "in") + unit(3*0.5, "mm"),
         y = 0.655,
         fontsize = 5, fontfamily = "Helvetica")

dev.off()

# CHROMOSOME BARPLOTS -----------------------------------------------------

ctl_sig_genes <- read_csv("data/sex_de/ctl_sexDE_pval01.csv") %>%
  mutate(condition = "PBS") %>%
  mutate(sex = ifelse(log2FoldChange < 0, "Female", "Male"))
fnf_sig_genes <- read_csv("data/sex_de/fnf_sexDE_pval01.csv") %>%
  mutate(condition = "FN-f") %>%
  mutate(sex = ifelse(log2FoldChange < 0, "Female", "Male"))

all_sig_sex_genes <- bind_rows(ctl_sig_genes,
                               fnf_sig_genes)
all_sig_sex_genes$seqnames <- factor(all_sig_sex_genes$seqnames,
                                     levels = c(as.character(1:22), "X", "Y"))
all_sig_sex_genes$condition <- factor(all_sig_sex_genes$condition,
                                      levels = c("PBS", "FN-f"))


ggplot(all_sig_sex_genes, aes(x = seqnames, fill = sex)) +
  geom_bar() +
  facet_wrap(vars(condition), ncol = 1, scales = "free") +
  scale_fill_manual(values = c(sexColors[["F"]], sexColors[["M"]])) +
  scale_x_discrete(name = "Chromosome", drop = FALSE) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 30),
                     labels = c(10, 20, 30),
                     breaks = c(10, 20, 30),
                     name = 'Number of 
                     <span style = "color:#4788BA;">**male**</span> and
                     <span style = "color:#DD8492;">**female**</span><br> sex-specific genes') +
  theme(legend.position = "none",
        text = element_text(family = "Helvetica"),
        panel.background = element_blank(),
        panel.spacing = unit(0, "cm"),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(linewidth = 0.25),
        axis.ticks.length.y = unit(-0.1, "cm"),
        axis.title.y = element_markdown(family = "Helvetica", size = 8),
        axis.title.x = element_text(size = 8),
        axis.line = element_line(linewidth = 0.25),
        axis.text.x = element_text(color = "black", size = 6,
                                   face = c(rep("plain", 22), rep("bold", 2))),
        axis.text.y = element_text(color = "black", size = 8),
        strip.background = element_blank(),
        strip.text = element_text(size = 10, vjust = 0, face = "bold"))


ggsave(filename = "plots/sexDE_Fig2/sex_fnfpbs_numGenes_barplot.pdf",
       width = 4, height = 3.5, units = "in")



# VENN DIAGRAM OF PBS AND FNF SEX GENES -----------------------------------
ctl_sig_genes <- read_csv("data/sex_de/ctl_sexDE_pval01.csv")
fnf_sig_genes <- read_csv("data/sex_de/fnf_sexDE_pval01.csv")


pbs_fnf_sex_genes <- tibble(values = unique(c(ctl_sig_genes$gene_id, fnf_sig_genes$gene_id))) %>%
  mutate(PBS = values %in% ctl_sig_genes$gene_id,
         FNF = values %in% fnf_sig_genes$gene_id)


venn_diagram <- ggplot(pbs_fnf_sex_genes, aes(A = PBS, B = FNF)) +
  geom_venn(set_names = c("PBS", "FN-f"), 
            fill_color = c(log2fcColors[["-"]], log2fcColors[["+"]]), 
            stroke_color = NA, auto_scale = TRUE, show_percentage = FALSE,
            text_size = 6, set_name_size = 6) +
  coord_fixed()  +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank())
venn_font(venn_diagram, font = "Helvetica")
ggsave(filename = "plots/sexDE_Fig2/sex_pbs_fnf_venndiagram.pdf", width = 4, height = 4, units = "in")

#### Checking directions of effect for each set

# Overlap
overlap_sex_genes <- ctl_sig_genes[which(ctl_sig_genes$gene_id %in% fnf_sig_genes$gene_id), "gene_id"]
ctl_overlap <- ctl_sig_genes %>%
  filter(gene_id %in% overlap_sex_genes$gene_id) %>%
  dplyr::select(gene_id, symbol, log2FoldChange) %>%
  dplyr::rename(ctl_l2fc = log2FoldChange)

fnf_overlap <- fnf_sig_genes %>%
  filter(gene_id %in% overlap_sex_genes$gene_id) %>%
  dplyr::select(gene_id, symbol, log2FoldChange) %>%
  dplyr::rename(fnf_l2fc = log2FoldChange)

overlap_directions <- left_join(ctl_overlap,
                                fnf_overlap,
                                by = c("gene_id", "symbol")) %>%
  mutate(concordant_dir = ifelse(sign(ctl_l2fc) == sign(fnf_l2fc), TRUE, FALSE))

overlap_percent_concordant <- length(which(overlap_directions$concordant_dir))/nrow(overlap_directions)

# ctl only, looking up in all fnf
ctl_unique <- ctl_sig_genes[which(!ctl_sig_genes$gene_id %in% fnf_sig_genes$gene_id), "gene_id"]

ctl_ctl_unique <- ctl_sig_genes %>%
  filter(gene_id %in% ctl_unique$gene_id) %>%
  dplyr::select(gene_id, symbol, log2FoldChange) %>%
  dplyr::rename(ctl_l2fc = log2FoldChange)

fnf_ctl_unique <- read_csv("data/sex_de/fnf_sex_shrink.csv",
         col_select = c("gene_id", "symbol", "log2FoldChange")) %>%
  filter(gene_id %in% ctl_unique$gene_id) %>%
  dplyr::rename(fnf_l2fc = log2FoldChange)
  
ctl_unique_directions <- left_join(ctl_ctl_unique,
                                   fnf_ctl_unique,
                                   by = c("gene_id", "symbol")) %>%
  mutate(concordant_dir = ifelse(sign(ctl_l2fc) == sign(fnf_l2fc), TRUE, FALSE))

ctl_unique_percent_concordant <- length(which(ctl_unique_directions$concordant_dir))/nrow(ctl_unique_directions)

# fnf only, looking up in all ctl
fnf_unique <- fnf_sig_genes[which(!fnf_sig_genes$gene_id %in% ctl_sig_genes$gene_id), "gene_id"]

fnf_fnf_unique <- fnf_sig_genes %>%
  filter(gene_id %in% fnf_unique$gene_id) %>%
  dplyr::select(gene_id, symbol, log2FoldChange) %>%
  dplyr::rename(fnf_l2fc = log2FoldChange)

ctl_fnf_unique <- read_csv("data/sex_de/ctl_sex_shrink.csv",
                           col_select = c("gene_id", "symbol", "log2FoldChange")) %>%
  filter(gene_id %in% fnf_unique$gene_id) %>%
  dplyr::rename(ctl_l2fc = log2FoldChange)

fnf_unique_directions <- left_join(fnf_fnf_unique,
                                   ctl_fnf_unique,
                                   by = c("gene_id", "symbol")) %>%
  mutate(concordant_dir = ifelse(sign(ctl_l2fc) == sign(fnf_l2fc), TRUE, FALSE))

fnf_unique_percent_concordant <- length(which(fnf_unique_directions$concordant_dir))/nrow(fnf_unique_directions)
# OVERLAP WITH FNF AND OA DIFF GENES --------------------------------------

# Union of sex degenes significant in control and FN-f
ctl_sex_degenes <- read_csv("data/sex_de/ctl_sexDE_pval01.csv")
fnf_sex_degenes <- read_csv("data/sex_de/fnf_sexDE_pval01.csv")

sex_degenes_union <- full_join(ctl_sex_degenes %>%
                                 dplyr::select("symbol", "gene_id", "log2FoldChange") %>%
                                 mutate(condition = "ctl"),
                               fnf_sex_degenes %>%
                                 dplyr::select("symbol", "gene_id", "log2FoldChange") %>%
                                 mutate(condition = "fnf")) %>%
  distinct(gene_id, .keep_all = TRUE) %>%
  dplyr::rename(sex_log2FoldChange = log2FoldChange)


# Differential FN-f, OA, and overlapping FN-f/OA genes
fnf_de_genes <- read_csv("data/condition_de/sig_deGenes_pval01_l2fc2.csv")
oa_de_genes <- read_csv("data/RAAK/RAAK_genes.csv",
                        col_select = c("ENSEMBL", "HGNC", "RAAK_PVAL",
                                       "RAAK_FC", "RAAK_LFC")) %>%
  dplyr::rename(symbol = HGNC) %>%
  dplyr::rename(log2FoldChange = RAAK_LFC)
fnf_oa_de_genes <- fnf_de_genes %>%
  filter(symbol %in% oa_de_genes$symbol)


# Get sex de genes only differential in FN-f
sex_fnf_degenes <- sex_degenes_union %>%
  filter(symbol %in% fnf_de_genes$symbol & !symbol %in% oa_de_genes$symbol) 

# All sex de genes differential in FN-f
sex_fnf_degenes_all <- sex_degenes_union %>%
  filter(symbol %in% fnf_de_genes$symbol) 

# Get sex de genes only differential in OA
sex_oa_degenes <- sex_degenes_union %>%
  filter(!symbol %in% fnf_de_genes$symbol & symbol %in% oa_de_genes$symbol) 

# All sex de genes differential in OA
sex_oa_degenes_all <- sex_degenes_union %>%
  filter(symbol %in% oa_de_genes$symbol) %>%
  left_join(oa_de_genes %>% dplyr::select(symbol, log2FoldChange)) %>%
  dplyr::rename(oa_log2FoldChange = log2FoldChange) %>%
  mutate(oa_group = ifelse(oa_log2FoldChange < 0, "down", "up"))

# Get sex de genes differential in FN-f and OA
sex_fnfoa_degenes <- sex_degenes_union %>%
  filter(symbol %in% fnf_de_genes$symbol & symbol %in% oa_de_genes$symbol) 

## Permutation tests for enrichment
de_genes_results <- read_csv("data/condition_de/de_genes_results.csv")

fnf_pval <- gene_enrichment_permTest(group = "fnf")
oa_pval <- gene_enrichment_permTest(group = "oa")
fnf_oa_pval <- gene_enrichment_permTest(group = "fnf_oa")

###### Sex-specific OA gene count boxplots in PBS and FNF
load("data/sex_de/dds_sex_ctl.rda")
load("data/sex_de/dds_sex_fnf.rda")

all_gene_sex_counts <- list()
for (geneRow in 1:nrow(sex_oa_degenes_all)){
  gene <- sex_oa_degenes_all[geneRow, ]
  
  ctl_gene_counts <- get_gene_sex_Counts(gene = gene[["gene_id"]], dds = dds_sex_ctl) %>%
    mutate(condition = "PBS")
  fnf_gene_counts <- get_gene_sex_Counts(gene = gene[["gene_id"]], dds = dds_sex_fnf) %>%
    mutate(condition = "FN-f")
  
  all_gene_counts <- bind_rows(ctl_gene_counts, fnf_gene_counts)
  all_gene_sex_counts[[gene[["gene_id"]]]] <- all_gene_counts
}

all_sex_oa_counts <- left_join(bind_rows(all_gene_sex_counts), 
                               sex_oa_degenes_all %>%
                                 dplyr::select(symbol, gene_id, oa_group), by = "gene_id")
all_sex_oa_counts$oa_group <- factor(all_sex_oa_counts$oa_group,
                                     levels = c("up", "down"))
all_sex_oa_counts$symbol <- factor(all_sex_oa_counts$symbol,
                                   levels = c("SERPINE2", "TNFRSF12A", "PDPN", "RARRES2"))
all_sex_oa_counts$condition <- factor(all_sex_oa_counts$condition,
                                      levels = c("PBS", "FN-f"))

gene_maxCounts <- all_sex_oa_counts %>%
  group_by(symbol) %>%
  summarise(maxCount = max(log2(count)))

ctl_sex_OAgenes <- read_csv("data/sex_de/ctl_sexDE_pval01.csv", 
                            col_select = c("gene_id", "symbol", "padj", "log2FoldChange")) %>%
  filter(gene_id %in% sex_oa_degenes_all$gene_id) %>%
  mutate(effect_dir = ifelse(log2FoldChange < 0, "F", "M")) %>%
  dplyr::select(-log2FoldChange) %>%
  dplyr::rename(PBS = padj)
fnf_sex_OAgenes <- read_csv("data/sex_de/fnf_sexDE_pval01.csv", 
                            col_select = c("gene_id", "symbol", "padj", "log2FoldChange")) %>%
  filter(gene_id %in% sex_oa_degenes_all$gene_id) %>%
  mutate(effect_dir = ifelse(log2FoldChange < 0, "F", "M")) %>%
  dplyr::select(-log2FoldChange) %>%
  dplyr::rename(`FN-f` = padj)

all_sig_sex_OAgenes <- full_join(ctl_sex_OAgenes, fnf_sex_OAgenes) %>%
  pivot_longer(cols = c("PBS", "FN-f"), names_to = "condition", values_to = "padj") %>%
  left_join(gene_maxCounts) %>%
  filter(!is.na(padj)) %>%
  mutate(signif = "*")
all_sig_sex_OAgenes <- merge(all_sig_sex_OAgenes, data.frame("Sex" = c("M", "F")))
all_sig_sex_OAgenes$condition <- factor(all_sig_sex_OAgenes$condition,
                                        levels = c("PBS", "FN-f"))

up_boxplots <- ggplot(all_sex_oa_counts %>% 
                        filter(oa_group == "up"), 
                      aes(x = Sex, y = log2(count), fill = Sex)) +
  geom_boxplot(outlier.shape = NA) +
  ggh4x::facet_grid2(condition ~ symbol, switch = "both", scales = "free",
                     independent = "all") +
  geom_line(data = all_sig_sex_OAgenes %>% 
              filter(symbol %in% c("SERPINE2", "TNFRSF12A")), 
            aes(x = Sex, y = maxCount + 0.5, group = symbol), inherit.aes = FALSE,
            linewidth = 0.5) +
  geom_star(data = all_sig_sex_OAgenes %>% 
              filter(symbol %in% c("SERPINE2", "TNFRSF12A")),
            aes(x = 1.5, y = Inf, starshape = signif, fill = effect_dir), inherit.aes = FALSE,
            size = 3) +
  scale_fill_manual(values = c(sexColors[["F"]], sexColors[["M"]])) +
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

down_boxplots <- ggplot(all_sex_oa_counts %>% 
                          filter(oa_group == "down"), 
                        aes(x = Sex, y = log2(count), fill = Sex)) +
  geom_boxplot(outlier.shape = NA) +
  ggh4x::facet_grid2(condition ~ symbol, switch = "both", scales = "free",
                     independent = "all") +
  geom_line(data = all_sig_sex_OAgenes %>% 
              filter(symbol %in% c("PDPN", "RARRES2")), 
            aes(x = Sex, y = maxCount + 0.5, group = symbol), inherit.aes = FALSE,
            linewidth = 0.5) +
  geom_star(data = all_sig_sex_OAgenes %>% 
              filter(symbol %in% c("PDPN", "RARRES2")),
            aes(x = 1.5, y = Inf, starshape = signif, fill = effect_dir), 
            inherit.aes = FALSE,
            size = 3) +
  scale_fill_manual(values = c(sexColors[["F"]], sexColors[["M"]])) +
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
        strip.text.y.left = element_blank(),
        strip.text.x.bottom = element_text(face = "bold", size = 12),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none",
        panel.spacing.y = unit(0.5, "cm"),
        plot.title = element_text(hjust = 0.5)) +
  ggtitle("Down in OA")

# Combine with patchwork
up_boxplots + down_boxplots


ggsave(filename = "plots/sexDE_Fig2/oa_sex_genes_boxplots.pdf", width = 9, height = 5, units= "in")


# Overlap with sex-biased genes from GTEx ---------------------------------

# Significant from PBS and FNF
ctl_sig_genes <- read_csv("data/sex_de/ctl_sexDE_pval01.csv", 
                          col_select = c("gene_id", "symbol", 
                                         "log2FoldChange", "lfcSE"))
fnf_sig_genes <- read_csv("data/sex_de/fnf_sexDE_pval01.csv",
                          col_select = c("gene_id", "symbol", 
                                         "log2FoldChange", "lfcSE"))
union_sig_genes <- union(ctl_sig_genes$gene_id, fnf_sig_genes$gene_id)

sex_genes <- bind_rows(ctl_sig_genes %>%
                         filter(gene_id %in% union_sig_genes),
                       fnf_sig_genes %>%
                         filter(gene_id %in% union_sig_genes)) %>%
  distinct(gene_id, .keep_all = TRUE) %>%
  mutate(sex = ifelse(log2FoldChange < 0, "female", "male"))

# Pull union list of our data
chond_sbgenes <- sex_genes %>%
  mutate(sex = ifelse(log2FoldChange < 0, "female", "male")) %>%
  mutate(tissue = "Chondrocytes") %>%
  dplyr::rename(effsize = log2FoldChange) %>%
  dplyr::rename(effsize_se = lfcSE)

# Read in gtex sex-biased genes
gtex_signif_sbgenes <- 
  read_delim("data/GTEx/GTEx_Analysis_v8_sbgenes/signif.sbgenes.txt") %>%
  # convert IDs to ones compatible with ours
  mutate(gene_id = gsub("\\..*", "", gene)) %>%
  # Flip sign of effsize to match ours (gtex positive is female and negative is male)
  mutate(effsize = -1*effsize) %>%
  mutate(sex = ifelse(effsize < 0, "female", "male")) %>%
  # filter for genes in chond_sbgenes
  filter(gene_id %in% chond_sbgenes$gene_id) %>%
  mutate(tissue = gsub("_", " ", tissue))


# Join with our data
all_data_sbgenes <- bind_rows(chond_sbgenes %>% 
                                dplyr::select(gene_id, tissue, sex, effsize, effsize_se),
                              gtex_signif_sbgenes %>% 
                                dplyr::select(gene_id, tissue, sex, effsize, effsize_se)) %>%
  complete(gene_id, tissue) %>% 
  left_join(chond_sbgenes %>% dplyr::select(gene_id, symbol), by = "gene_id")

# Order genes by how many overlap
geneOverlaps <- all_data_sbgenes %>%
  group_by(symbol) %>%
  filter(!is.na(effsize)) %>%
  summarize(nOverlap = dplyr::n()) %>%
  arrange(nOverlap)

all_data_sbgenes <- left_join(all_data_sbgenes, geneOverlaps, by = "symbol") %>%
  group_by(nOverlap) %>%
  arrange(sex, .by_group = TRUE) %>%
  ungroup()

all_data_sbgenes$tissue <- factor(all_data_sbgenes$tissue,
                                  levels = c("Chondrocytes", 
                                             unique(gtex_signif_sbgenes$tissue)))
all_data_sbgenes$symbol <- factor(all_data_sbgenes$symbol,
                                  levels = unique(all_data_sbgenes$symbol))

# Add fake data to create alternating colored background for rows
all_data_sbgenes <- all_data_sbgenes %>%
  mutate(xmin = 0.5, xmax = 45.5, 
         y_position = as.numeric(symbol),
         ymin = y_position - 0.5,
         ymax = y_position + 0.5) %>%
  pivot_longer(cols = c(xmin, xmax), values_to = "x", names_to = "xmin_xmax") %>%
  dplyr::select(-xmin_xmax) %>%
  mutate(fill = ifelse(y_position %% 2 == 0, "a", "b"))

# Create fontFaces vector for bolding chond specific genes
fontFaces <- all_data_sbgenes %>%
  reframe(face = ifelse(nOverlap == 1, "bold", "plain"), .by = symbol) %>%
  distinct()


# Get log2FC for all sig sex genes in PBS and FNF
ctl_sex_results <- read_csv("data/sex_de/ctl_sex_shrink.csv", 
                            col_select = c("gene_id", "symbol", "log2FoldChange", "lfcSE")) %>%
  filter(gene_id %in% sex_genes$gene_id) %>%
  mutate(condition = "PBS") %>%
  mutate(sex = ifelse(log2FoldChange < 0, "female", "male"))
fnf_sex_results <- read_csv("data/sex_de/fnf_sex_shrink.csv",
                            col_select = c("gene_id", "symbol", "log2FoldChange", "lfcSE")) %>%
  filter(gene_id %in% sex_genes$gene_id) %>%
  mutate(condition = "FN-f") %>%
  mutate(sex = ifelse(log2FoldChange < 0, "female", "male"))
all_sex_results <- bind_rows(ctl_sex_results,
                             fnf_sex_results) %>%
  complete(symbol, condition)
# Put genes in same order as heatmap
all_sex_results$symbol <- factor(all_sex_results$symbol, 
                                 levels = unique(all_data_sbgenes$symbol))
all_sex_results$condition <- factor(all_sex_results$condition,
                                    levels = c("PBS", "FN-f"))

# Add fake data for stripes
all_sex_results <- all_sex_results %>%
  mutate(xmin = -5, xmax = 20, 
         y_position = as.numeric(symbol),
         ymin = y_position - 0.5,
         ymax = y_position + 0.5) %>%
  pivot_longer(cols = c(xmin, xmax), values_to = "x", names_to = "xmin_xmax") %>%
  dplyr::select(-xmin_xmax) %>%
  mutate(fill = ifelse(y_position %% 2 == 0, "a", "b"))


# Median l2fc across GTEx tissues
gtex_median_l2fc <- right_join(gtex_signif_sbgenes, 
                            chond_sbgenes %>% dplyr::select(gene_id, symbol), 
                            by = "gene_id") %>%
  group_by(symbol) %>%
  summarise(median_l2fc = median(effsize), 
            q1 = quantile(effsize, probs = 0.25, na.rm = TRUE),
            q3 = quantile(effsize, probs = 0.75, na.rm = TRUE)) %>%
  mutate(sex = ifelse(median_l2fc < 0, "female", "male"))
# Put genes in same order as heatmap
gtex_median_l2fc$symbol <- factor(gtex_median_l2fc$symbol,
                                levels = unique(all_data_sbgenes$symbol))

# Add fake data for stripes
gtex_median_l2fc <- gtex_median_l2fc %>%
  mutate(xmin = -5, xmax = 20, 
         y_position = as.numeric(symbol),
         ymin = y_position - 0.5,
         ymax = y_position + 0.5) %>%
  pivot_longer(cols = c(xmin, xmax), values_to = "x", names_to = "xmin_xmax") %>%
  dplyr::select(-xmin_xmax) %>%
  mutate(fill = ifelse(y_position %% 2 == 0, "a", "b"))


##### Heatmap of M vs F for each gene across tissues
sex_cell_heatmap <- ggplot(all_data_sbgenes, aes(x = tissue, y = y_position)) +
  scale_x_discrete(position = "top") +
  coord_cartesian(clip = "off") +
  geom_ribbon(aes(x = x, ymin = ymin, ymax = ymax, group = y_position, fill = fill), 
              inherit.aes = FALSE, alpha = 0.2) +
  scale_fill_manual(values = c("grey50", "grey75")) +
  ggnewscale::new_scale_fill() +
  geom_tile(aes(fill = sex)) +
  scale_fill_manual(values = c(sexColors[["F"]], sexColors[["M"]]), na.value = NA) +
  scale_y_continuous(breaks = unique(all_data_sbgenes$y_position), 
                     labels = unique(all_data_sbgenes$symbol), expand = c(0, 0)) +
  theme(panel.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 0, color = "black", size = 10,
                                   face = c("bold", rep("plain", 44))),
        axis.text.y = element_text(color = "black", vjust = 0.5,
                                   face = fontFaces$face),
        axis.title = element_blank(),
        text = element_text(family = "Helvetica"),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(color = "grey25", linewidth = 0.2),
        axis.ticks.length.y = unit(2, "mm"),
        legend.position = "none",
        plot.margin = margin(10, 50, 10, 10))


### log2FC of these genes in chondrocytes (PBS and FNF)
cellOverlap_chond_l2fc <- ggplot(all_sex_results, aes(x = log2FoldChange, y = y_position)) +
  geom_ribbon(aes(x = x, ymin = ymin, ymax = ymax, group = y_position, fill = fill), 
              inherit.aes = FALSE, alpha = 0.2) +
  scale_fill_manual(values = c("grey50", "grey75")) +
  geom_vline(xintercept = 0, lty = 2) +
  geom_segment(aes(x = log2FoldChange-1.96*lfcSE, xend = log2FoldChange+1.96*lfcSE,
                   y = y_position, yend = y_position), color = "grey25") +
  geom_point(aes(color = sex)) +
  # Add NA text for missing points +
  geom_text(data = all_sex_results %>%
              filter(is.na(log2FoldChange)) %>%
              mutate(log2FoldChange = "NA"),
            aes(x = 1, y = y_position, label = log2FoldChange), 
            inherit.aes = FALSE, size = 3, hjust = 0, family = "Helvetica") +
  facet_wrap(vars(condition)) +
  coord_cartesian(clip = "off") +
  scale_color_manual(values = c(sexColors[["F"]], sexColors[["M"]])) +
  scale_x_continuous(limits = c(-5, 20), name = "Gene sex log~2~(fold change)<br> in chondrocytes",
                     expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0))+
  theme(axis.title.x = element_markdown(family = "Helvetica", size = 9),
        axis.title.y = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 10),
        legend.position = "none",
        panel.spacing.x = unit(1, "cm"),
        panel.border = element_rect(fill = NA, color = "grey25"),
        panel.background = element_blank(),
        axis.text.x = element_text(color = "black", size = 10),
        axis.text.y = element_blank(),
        axis.line.x = element_line(color = "grey25"),
        axis.ticks = element_blank(),
        text = element_text(family = "Helvetica"),
        plot.margin = margin(10, 10, 10, 10))

#### Median log2FC of genes across GTEx tissues

cellOverlap_gtex_l2fc <- ggplot(gtex_median_l2fc, aes(x = median_l2fc, y = y_position)) +
  geom_ribbon(aes(x = x, ymin = ymin, ymax = ymax, group = y_position, fill = fill),
              inherit.aes = FALSE, alpha = 0.2) +
  scale_fill_manual(values = c("grey50", "grey75")) +
  geom_vline(xintercept = 0, lty = 2) +
  geom_segment(aes(x = q1, xend = q3,
                   y = y_position, yend = y_position), color = "grey25") +
  geom_point(aes(color = sex)) +
  coord_cartesian(clip = "off") +
  scale_color_manual(values = c(sexColors[["F"]], sexColors[["M"]])) +
  scale_x_continuous(expand = c(0, 0), name = "Median log~2~(fold change)<br>across GTEx tissues") +
  scale_y_continuous(expand = c(0, 0)) +
  theme(axis.title.x = element_markdown(family = "Helvetica", size = 9),
        axis.title.y = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 10),
        legend.position = "none",
        panel.spacing.x = unit(1, "cm"),
        panel.border = element_rect(fill = NA, color = "grey25"),
        panel.background = element_blank(),
        axis.text.x = element_text(color = "black", size = 10),
        axis.text.y = element_blank(),
        axis.line.x = element_line(color = "grey25"),
        axis.ticks = element_blank(),
        text = element_text(family = "Helvetica"),
        plot.margin = margin(10, 20, 10, 20))

sex_cell_heatmap  + cellOverlap_gtex_l2fc + cellOverlap_chond_l2fc +
  plot_layout(ncol = 3, widths = c(4, 0.55, 1.28))

ggsave(filename = "plots/sexDE_Fig2/sex_celltypeOverlap_heatmap_gtexl2fc_chondl2fc.pdf",
       width = 15, height = 15, units = "in", bg = "transparent")
