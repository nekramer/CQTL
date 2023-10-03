library(ComplexHeatmap)
library(DESeq2)
library(tidyverse)
library(circlize)
library(RColorBrewer)
library(plotgardener)
library(grid)
library(rrvgo)
library(ggnewscale)
library(patchwork)
library(rsvg)
library(grImport2)
library(org.Hs.eg.db)
source("../plotting_utils.R")


# Functions ---------------------------------------------------------------

get_sample_l2fc <- function(gene, countMatrix){
  
  # Extract row of gene from countMatrix
  gene_counts <- countMatrix[gene,]
  
  # Convert to dataframe and extract donors/conditions into separate columns
  donor_gene_counts <- data.frame(gene_counts) %>%
    rownames_to_column(var = "Sample") %>%
    separate_wider_delim("Sample", 
                         delim = "_", 
                         names = c(NA, "Donor", NA, "Condition", NA, NA)) %>%
    # Group by each donor and calculate l2FC 
    group_by(Donor) %>%
    summarize(log2FC = 
                log2(gene_counts[Condition == "FNF"]/gene_counts[Condition == "CTL"])) %>%
    ungroup() %>%
    mutate(ENSEMBL = gene)
  
  return(donor_gene_counts)
}

# HEATMAP -----------------------------------------------------------------

load("data/condition_de/differential_expression_dds.rda")
sig_degenes <- 
  read_csv("data/condition_de/sig_deGenes_pval01_l2fc2.csv") %>%
  mutate(log2FC_dir = ifelse(log2FoldChange < 0, "-", "+")) %>%
  arrange(log2FC_dir)
donorSamplesheet <- read_csv("data/donorSamplesheet.csv")

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
annotations$Age <- cut(annotations$Age, breaks = seq(30, 90, 10),
                       labels = c("31-40", "41-50", "51-60", "61-70", "71-80", "81-90"))
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
  col = list(Age = c("31-40" = ageColors[6],
                     "41-50" = ageColors[5],
                     "51-60" = ageColors[4],
                     "61-70" = ageColors[3],
                     "71-80" = ageColors[2],
                     "81-90" = ageColors[1]),
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
new_col_order <- c(col_order[103:204], col_order[1:102])

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
pdf(file = "plots/conditionDE_Fig1/DE_heatmap.pdf", width = 6.5, height = 5.25)
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
plotText(label = "31", 
         x = 5.75,
         y = 0.18,
         fontfamily = "Helvetica", fontsize = 4, just = c("left", "top"))
plotText(label = "41", 
         x = unit(5.75, "in") + unit(3, "mm"),
         y = 0.18,
         fontfamily = "Helvetica", fontsize = 4, just = c("left", "top"))
plotText(label = "51", 
         x = unit(5.75, "in") + unit(3*2, "mm"),
         y = 0.18,
         fontfamily = "Helvetica", fontsize = 4, just = c("left", "top"))
plotText(label = "61", 
         x = unit(5.75, "in") + unit(3*3, "mm"),
         y = 0.18,
         fontfamily = "Helvetica", fontsize = 4, just = c("left", "top"))
plotText(label = "71", 
         x = unit(5.75, "in") + unit(3*4, "mm"),
         y = 0.18,
         fontfamily = "Helvetica", fontsize = 4, just = c("left", "top"))
plotText(label = "81", 
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

# GO AND KEGG BARPLOTS ----------------------------------------------------
#### GO

# Get reduced, significant GO terms for each category
upsig_go_data <- 
  read_delim("data/homer/homer_upsig_deGenes_pval01_l2fc2/biological_process.txt") %>%
  mutate(pval = exp(1)^logP) %>%
  filter(pval < 0.01)
upsig_go <- reduceGO(upsig_go_data,
                     category = "Upregulated")


downsig_go_data <- 
  read_delim("data/homer/homer_downsig_deGenes_pval01_l2fc2/biological_process.txt") %>%
  mutate(pval = exp(1)^logP) %>%
  filter(pval < 0.01)
downsig_go <- reduceGO(downsig_go_data,
                       category = "Downregulated")


# Select 5 each for plotting
upsig_go_plotting <- upsig_go %>%
  filter(parentTerm %in% c("response to cytokine", "cell surface receptor signaling pathway", "collagen catabolic process", "regulation of cell-cell adhesion", "immune effector process")) |> 
  arrange(`-log10pval`)

downsig_go_plotting <- downsig_go %>%
  filter(parentTerm %in% c("anatomical structure development", 
                           "skeletal system development",
                           "regulation of developmental process", 
                           "positive regulation of collagen biosynthetic process",
                           "regulation of multicellular organismal process")) %>%
  arrange(`-log10pval`)

# Combine into one
go_plotting <- bind_rows(upsig_go_plotting, downsig_go_plotting)
go_plotting$parentTerm <- factor(go_plotting$parentTerm, levels = go_plotting$parentTerm) 
go_plotting$category <- factor(go_plotting$category, levels = c("Upregulated", "Downregulated"))

# Plot all in barplot
ggplot(go_plotting, aes(x = `-log10pval`, y = parentTerm, fill = category)) +
  geom_vline(xintercept = 10, color = "grey75", alpha = 0.4) +
  geom_vline(xintercept = 20, color = "grey75", alpha = 0.4) +
  geom_vline(xintercept = 30, color = "grey75", alpha = 0.4) +
  geom_vline(xintercept = 40, color = "grey75", alpha = 0.4) +
  geom_vline(xintercept = 50, color = "grey75", alpha = 0.4) +
  geom_bar(stat = "identity") +
  scale_x_continuous(limits = c(0, 51), expand = c(0, 0), name = "-log~10~pval",
                     breaks = seq(0, 50, 10)) +
  scale_fill_manual(values = c(log2fcColors[["+"]], log2fcColors[["-"]])) +
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
        axis.ticks = element_blank(),
        axis.line.x = element_line(),
        strip.text = element_text(size = 14, color = "black"),
        panel.spacing = unit(0, "mm"), 
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14)) +
  ggtitle("GO Terms")

ggsave(filename = "plots/conditionDE_Fig1/GO_barplots.pdf", width = 5, height = 8, units = "in")


#### KEGG

# Plot top 20 significant for each category
upsig_kegg_data <- read_delim("data/homer/homer_upsig_deGenes_pval01_l2fc2/kegg.txt") %>%
  mutate(pval = exp(1)^logP) %>%
  filter(pval < 0.01) %>%
  distinct(Term, .keep_all = TRUE) %>%
  mutate(`-log10pval` = -log10(pval)) %>%
  mutate(category = "Upregulated")
downsig_kegg_data <- read_delim("data/homer/homer_downsig_deGenes_pval01_l2fc2/kegg.txt") %>%
  mutate(pval = exp(1)^logP) %>%
  filter(pval < 0.01) %>%
  distinct(Term, .keep_all = TRUE) %>%
  mutate(`-log10pval` = -log10(pval)) %>%
  mutate(category = "Downregulated")


upsig_kegg_plotting <- upsig_kegg_data %>%
  filter(Term %in% c("TNF signaling pathway", "IL-17 signaling pathway", 
                     "Cytokine-cytokine receptor interaction", 
                     "NF-kappa B signaling pathway", "NOD-like receptor signaling pathway")) %>%
  arrange(`-log10pval`)

down_kegg_plotting <- downsig_kegg_data %>%
  filter(Term %in% c("Rap1 signaling pathway", "Drug metabolism - cytochrome P450",
                     "Calcium signaling pathway", "PPAR signaling pathway",
                     "cAMP signaling pathway")) %>%
  arrange(`-log10pval`)

kegg_plotting <- bind_rows(upsig_kegg_plotting, down_kegg_plotting)
kegg_plotting$Term <- factor(kegg_plotting$Term, levels = kegg_plotting$Term) 
kegg_plotting$category <- factor(kegg_plotting$category, levels = c("Upregulated", "Downregulated"))

ggplot(kegg_plotting, aes(x = `-log10pval`, y = Term, fill = category)) +
  geom_vline(xintercept = 5, color = "grey75", alpha = 0.4) +
  geom_vline(xintercept = 10, color = "grey75", alpha = 0.4) +
  geom_vline(xintercept = 15, color = "grey75", alpha = 0.4) +
  geom_vline(xintercept = 20, color = "grey75", alpha = 0.4) +
  geom_vline(xintercept = 25, color = "grey75", alpha = 0.4) +
  geom_bar(stat = "identity") +
  scale_x_continuous(expand = c(0, 0), name = "-log~10~pval", limits = c(0, 26),
                     breaks = seq(0, 25, 5)) +
  scale_fill_manual(values = c(log2fcColors[["+"]], log2fcColors[["-"]])) +
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

ggsave(filename = "plots/conditionDE_Fig1/KEGG_barplots.pdf", width = 5, height = 8, units = "in")

# TF MOTIFS AND TF GENE EXPRESSION ----------------------------------------

de_genes_results <- read_csv("data/condition_de/de_genes_results.csv",
                             col_select = c("symbol", "log2FoldChange"))
load("data/condition_de/differential_expression_dds.rda")
normCounts <- counts(dds, normalized = TRUE)

### Upregulated motifs

# Read in and subset
upsig_knownmotifs <- read_delim("data/homer/homer_upsig_deGenes_pval01_l2fc2/knownResults.txt") %>%
  # Convert percentages to numbers
  mutate(across(c(`% of Target Sequences with Motif`, 
                  `% of Background Sequences with Motif`),
                ~ gsub("%", "", .))) %>%
  mutate(across(c(`% of Target Sequences with Motif`, 
                  `% of Background Sequences with Motif`), as.numeric)) %>%
  # Calculate log2 enrichment
  mutate(log2enrichment = log2(`% of Target Sequences with Motif`/`% of Background Sequences with Motif`)) %>%
  # Calculate -log10pval
  mutate(log10pval = -log10(exp(`Log P-value`))) %>%
  slice_max(order_by = log10pval, n = 4) %>%
  mutate(motifLogo = paste0("data/homer/homer_upsig_deGenes_pval01_l2fc2/knownResults/known", row_number(), ".logo.svg"))


# Pull out first part of motif name
upsig_knownmotifs$Name <- unlist(lapply(str_split(upsig_knownmotifs$`Motif Name`, 
                                                  '[(]'), `[[`, 1))

# Remove second NFKB
upsig_knownmotifs <- upsig_knownmotifs[-2,]

# Based on motifs, assign genes (and family members) that encode the TF
genes <- bind_rows(expand_grid(Name = upsig_knownmotifs$Name[1], 
                               Gene = c("NFKB1", "NFKB2", "RELA", "RELB", "REL")),
                   expand_grid(Name = upsig_knownmotifs$Name[2],
                               Gene = c("FOSL2", "FOS", "FOSB", "FOSL1")),
                   expand_grid(Name = upsig_knownmotifs$Name[3],
                               Gene = c("JUN", "JUNB", "JUND")))

upsig_knownmotifs <- left_join(upsig_knownmotifs, 
                               genes, 
                               by = "Name", 
                               multiple = "all")

# Get ENSEMBL IDs of genes
gene_symbols <- AnnotationDbi::select(org.Hs.eg.db, 
                                      keys = upsig_knownmotifs$Gene, 
                                      keytype = "SYMBOL", 
                                      columns = c("SYMBOL", "ENSEMBL"))

upsig_knownmotifs <- left_join(upsig_knownmotifs,
                               gene_symbols,
                               by = join_by(Gene == SYMBOL))

# For each gene, calculate each donor's FNF/CTL log2FC
upsig_donor_l2fcs <- bind_rows(lapply(unique(upsig_knownmotifs$ENSEMBL), 
                                      get_sample_l2fc, 
                                      countMatrix = normCounts))

# Join with motif data
upsig_knownmotifs_l2fc  <- left_join(upsig_knownmotifs, 
                                     upsig_donor_l2fcs, by = "ENSEMBL") %>%
  mutate(category = "up")

### Downregulated motifs

# Read in and subset for top 3
downsig_knownmotifs <- read_delim("data/homer/homer_downsig_deGenes_pval01_l2fc2/knownResults.txt") %>%
  # Convert percentages to numbers
  mutate(across(c(`% of Target Sequences with Motif`, 
                  `% of Background Sequences with Motif`),
                ~ gsub("%", "", .))) %>%
  mutate(across(c(`% of Target Sequences with Motif`, 
                  `% of Background Sequences with Motif`), as.numeric)) %>%
  # Calculate log2 enrichment
  mutate(log2enrichment = log2(`% of Target Sequences with Motif`/`% of Background Sequences with Motif`)) %>%
  # Calculate -log10pval
  mutate(log10pval = -log10(exp(`Log P-value`))) %>%
  slice_max(order_by = log10pval, n = 3) %>%
  mutate(motifLogo = paste0("data/homer/homer_downsig_deGenes_pval01_l2fc2/knownResults/known", row_number(), ".logo.svg"))

# Pull out first part of motif name
downsig_knownmotifs$Name <- unlist(lapply(str_split(downsig_knownmotifs$`Motif Name`, 
                                                    '[(]'), `[[`, 1))

# Assign gene names
downsig_knownmotifs <- downsig_knownmotifs %>% 
  filter(Name == "Mef2c") %>%
  mutate(Gene = toupper(Name))


# Get ENSEMBL IDs of genes
gene_symbols <- AnnotationDbi::select(org.Hs.eg.db, 
                                      keys = downsig_knownmotifs$Gene, 
                                      keytype = "SYMBOL", 
                                      columns = c("SYMBOL", "ENSEMBL"))

downsig_knownmotifs <- left_join(downsig_knownmotifs,
                                 gene_symbols,
                                 by = join_by(Gene == SYMBOL))

# For each gene, calculate each donor's FNF/CTL log2FC
downsig_donor_l2fcs <- bind_rows(lapply(unique(downsig_knownmotifs$ENSEMBL), 
                                        get_sample_l2fc, 
                                        countMatrix = normCounts))

# Join with motif data
downsig_knownmotifs_l2fc  <- left_join(downsig_knownmotifs, 
                                       downsig_donor_l2fcs, by = "ENSEMBL") %>%
  mutate(category = "down")

# Combine data into one 
sig_knownmotifs_l2fc <- bind_rows(upsig_knownmotifs_l2fc,
                                  downsig_knownmotifs_l2fc)
sig_knownmotifs_l2fc$Name <- factor(sig_knownmotifs_l2fc$Name, 
                                    levels = c("NFkB-p65-Rel", "Fos", "JunB", "Mef2c"))
sig_knownmotifs_l2fc$category <- factor(sig_knownmotifs_l2fc$category, levels = c("up", "down"))
sig_knownmotifs_l2fc$Gene <- factor(sig_knownmotifs_l2fc$Gene, 
                                    levels = c("MEF2C", "JUND", "JUNB",
                                               "JUN", "FOS", "FOSL2", "FOSB",
                                               "FOSL1", "REL", "RELA", "NFKB1",
                                               "RELB", "NFKB2"))
tf_gene_plot <- ggplot(sig_knownmotifs_l2fc, aes(x = log2FC, y = Gene, color = category)) +
  geom_vline(xintercept = 4, color = "grey90") +
  geom_vline(xintercept = 3, color = "grey90") +
  geom_vline(xintercept = 2, color = "grey90") +
  geom_vline(xintercept = 1, color = "grey90") +
  geom_vline(xintercept = -1, color = "grey90") +
  geom_vline(xintercept = -2, color = "grey90") +
  geom_vline(xintercept = -3, color = "grey90") +
  geom_vline(xintercept = -4, color = "grey90") +
  geom_violin(color = NA, fill = "grey25", alpha = 0.2) +
  geom_jitter(size = 0.75) +
  geom_violin(color = "grey25", fill = NA, linewidth = 0.25) +
  stat_summary(fun = "median", geom = "crossbar", width = 0.5, color = "grey25", linewidth = 0.5) +
  stat_boxplot(geom = "errorbar", width = 0.25, linewidth = 0.5, color = "grey25") +
  geom_hline(data = tibble(Name = factor(c("NFkB-p65-Rel"), 
                                         levels = c("NFkB-p65-Rel", "Fos", "JunB", "Mef2c") )), aes(yintercept = Inf)) +
  geom_hline(data = tibble(Name = factor(c("Mef2c", "JunB", "Fos"), 
                                         levels = c("NFkB-p65-Rel", "Fos", "JunB", "Mef2c") )), aes(yintercept = Inf),
             color = "grey") +
  geom_vline(data = tibble(Name = factor(c("NFkB-p65-Rel", "Fos", "JunB", "Mef2c")),
                           levels = c("NFkB-p65-Rel", "Fos", "JunB", "Mef2c")), aes(xintercept = Inf)) +
  geom_vline(xintercept = 0, lty = 2) +
  
  scale_color_manual(values = c(log2fcColors[["+"]], log2fcColors[["-"]])) +
  scale_x_continuous(limits = c(-5, 5), expand = c(0, 0), 
                     breaks = c(-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5),
                     name = "log~2~(fold change)") +
  coord_cartesian(clip = "off") +
  facet_wrap(~Name, ncol = 1, scales = "free_y") +
  theme(panel.background = element_blank(),
        text = element_text(family = "Helvetica"),
        legend.position = "None",
        axis.ticks = element_blank(),
        panel.spacing = unit(0, "mm"),
        axis.title.y = element_blank(),
        axis.title.x = element_markdown(),
        strip.text = element_blank(),
        axis.line = element_line(),
        axis.text = element_text(color = "black", size = 10))

# Convert to gtable 
tf_gene_plot_table <- ggplot_gtable(ggplot_build(tf_gene_plot))

# Resize panels to reflect number of genes in each category
panel1 <- tf_gene_plot_table$layout$t[grep('panel-1-1', tf_gene_plot_table$layout$name)]
panel2 <- tf_gene_plot_table$layout$t[grep('panel-1-2', tf_gene_plot_table$layout$name)]
panel3 <- tf_gene_plot_table$layout$t[grep('panel-1-3', tf_gene_plot_table$layout$name)]
panel4 <- tf_gene_plot_table$layout$t[grep('panel-1-4', tf_gene_plot_table$layout$name)]

tf_gene_plot_table$heights[panel1] <- 1.25 * tf_gene_plot_table$heights[panel1]
tf_gene_plot_table$heights[panel2] <- 1 * tf_gene_plot_table$heights[panel2]
tf_gene_plot_table$heights[panel3] <- 0.75 * tf_gene_plot_table$heights[panel3]
tf_gene_plot_table$heights[panel4] <- 0.25 * tf_gene_plot_table$heights[panel4]


motifImages <- list()
# Make motif images with name and pvalues
for (motif in unique(sig_knownmotifs_l2fc$Name)){
  
  # Get image
  motifImg <- pictureGrob(readPicture(rawToChar(rsvg::rsvg_svg(sig_knownmotifs_l2fc %>%
                                                                 filter(Name == motif) %>%
                                                                 pull(motifLogo) %>%
                                                                 unique()))))
  
  motifGrab <- grid.grabExpr(expr = {
    grid.newpage()
    grid.draw(motifImg)
  })
  
  motifImages[[motif]] <- motifGrab
  
}

grDevices::cairo_pdf(paste0("plots/conditionDE_Fig1/tfmotif_expression.pdf"), 
                     width = 6, height = 6)
pageCreate(width = 6, height = 6, showGuides = FALSE)
plotGG(tf_gene_plot_table, x = 2, y = 0, width = 4, height = 6)
plotGG(motifImages$`NFkB-p65-Rel`, x = 2, y = -0.6, width = 1.75, height = 1.75, just = c("right", "top"))

plotText("NF-\u03BAB", x = 1.1, y = 0.5, just = "top",
         fontfamily = "Helvetica", fontsize = 11, fontface = "bold")

plotText("Upregulated", x = 1.1, y = 0.65, just = "top",
         fontsize = 11, fontfamily = "Helvetica")

plotText(paste0("pval = ", sig_knownmotifs_l2fc %>%
                  filter(Name == "NFkB-p65-Rel") %>%
                  pull(`P-value`) %>%
                  unique()),
         x = 1.1, y = 0.8, just = "top", 
         fontsize = 11, fontfamily = "Helvetica")


plotGG(motifImages$Fos, x = 2, y = 1.5, width = 1.75, height = 1.75, just = c("right", "top"))
plotText("Fos", x = 1.1, y = 2.55, just = "top",
         fontfamily = "Helvetica", fontsize = 11, fontface = "bold")
plotText("Upregulated", x = 1.1, y = 2.7, just = "top",
         fontsize = 11, fontfamily = "Helvetica")
plotText(paste0("pval = ", sig_knownmotifs_l2fc %>%
                  filter(Name == "Fos") %>%
                  pull(`P-value`) %>%
                  unique()),
         x = 1.1, y = 2.85, just = "top",
         fontsize = 11, fontfamily = "Helvetica")

plotGG(motifImages$JunB, x = 2, y = 3.2, width = 1.75, height = 1.75, just = c("right", "top"))
plotText("JunB", x = 1.1, y = 4.3, just = "top",
         fontfamily = "Helvetica", fontsize = 11, fontface = "bold")
plotText("Upregulated", x = 1.1, y = 4.45, just = "top",
         fontsize = 11, fontfamily = "Helvetica")
plotText(paste0("pval = ", sig_knownmotifs_l2fc %>%
                  filter(Name == "JunB") %>%
                  pull(`P-value`) %>%
                  unique()),
         x = 1.1, y = 4.6, just = "top",
         fontsize = 11, fontfamily = "Helvetica")


plotGG(motifImages$Mef2c, x = 2, y = 4.3, width = 1.75, height = 1.75, just = c("right", "top"))
plotText("Mef2c", x = 1.1, y = 5.35, just = "top",
         fontfamily = "Helvetica", fontsize = 11, fontface = "bold")
plotText("Downregulated", x = 1.1, y = 5.5, just = "top",
         fontsize = 11, fontfamily = "Helvetica")
plotText(paste0("pval = ", sig_knownmotifs_l2fc %>%
                  filter(Name == "Mef2c") %>%
                  pull(`P-value`) %>%
                  unique()),
         x = 1.1, y = 5.65, just = "top",
         fontsize = 11, fontfamily = "Helvetica")
dev.off()

# COMPARISON WITH OA FROM RAAK STUDY --------------------------------------
raak_genes <- read_csv("data/RAAK/RAAK_genes.csv",
                       col_select = c("ENSEMBL", "HGNC", "RAAK_PVAL",
                                      "RAAK_FC", "RAAK_LFC"))

fnf_genes <- read_csv("data/condition_de/de_genes_results.csv")

raak_up <- raak_genes %>%
  filter(RAAK_LFC > 0)

raak_down <- raak_genes %>%
  filter(RAAK_LFC < 0)


fnf_up_raak_subset <- fnf_genes %>%
  filter(gene_id %in% raak_up$ENSEMBL | symbol %in% raak_up$HGNC) %>%
  mutate(group = "Up in OA")

fnf_down_raak_subset <- fnf_genes %>%
  filter(gene_id %in% raak_down$ENSEMBL | symbol %in% raak_down$HGNC) %>%
  mutate(group = "Down in OA")


up_test <- wilcox.test(x = fnf_up_raak_subset$log2FoldChange,
                       y = fnf_genes %>% 
                         filter(!gene_id %in% fnf_up_raak_subset$gene_id) %>% 
                         pull(log2FoldChange),
                       alternative = "greater")

down_test <- wilcox.test(x = fnf_down_raak_subset$log2FoldChange,
                         y = fnf_genes %>% 
                           filter(!gene_id %in% fnf_down_raak_subset$gene_id) %>% 
                           pull(log2FoldChange),
                         alternative = "less")

fnf_genes_raak_subset <- bind_rows(fnf_up_raak_subset,
                                   fnf_down_raak_subset)
fnf_genes_raak_subset$group <- factor(fnf_genes_raak_subset$group, 
                                      levels = c("Up in OA", "Down in OA"))

ggplot(fnf_genes_raak_subset, aes(x = group, y = log2FoldChange, 
                                  color = group, fill = group)) +
  scale_fill_manual(values = c(log2fcColors[["+"]], log2fcColors[["-"]])) +
  scale_color_manual(values = c(darken(log2fcColors[["+"]], 0.3), 
                                darken(log2fcColors[["-"]], 0.3))) +
  scale_y_continuous(name = "log~2~(fold change)<br>in response to FN-f",
                     limits = c(-8, 8),
                     breaks = seq(-8, 8, 2),
                     expand = c(0, 0)) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_boxplot(outlier.shape = 1, width = 0.35) +
  stat_boxplot(geom = "errorbar", width = 0.2) +
  annotate("text",
           label = paste0("pval = ", format(up_test$p.value, digits = 3)),
           x = 1, y = Inf, family = "Helvetica", size = 3.5, vjust = 1) +
  annotate("text",
           label = paste0("pval = ",format(down_test$p.value, digits = 3)),
           x = 2, y = Inf, family = "Helvetica", size = 3.5, vjust = 1) +
  theme(text = element_text(family = "Helvetica"),
        panel.background = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_markdown(size = 10),
        legend.position = "None",
        axis.line = element_line(),
        axis.text.x = element_text(color = "black", size = 10),
        axis.text.y = element_text(color = "black", size = 10),
        axis.ticks.y = element_line(color = "black"),
        axis.ticks.length.y = unit(-0.1, "cm"))

ggsave(file = "plots/conditionDE_Fig1/RAAKOA_FNF_boxplots.pdf", width = 6, height = 6, units = "in")
