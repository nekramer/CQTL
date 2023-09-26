library(tidyverse)
library(ggtext)
library(grid)
library(rsvg)
library(grImport2)
library(plotgardener)
library(org.Hs.eg.db)
source("scripts/plotting_utils.R")

de_genes_results <- read_csv("data/de_genes_results.csv",
                             col_select = c("symbol", "log2FoldChange"))
load("data/differential_expression_dds.rda")
normCounts <- counts(dds, normalized = TRUE)

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


# Upregulated motifs ------------------------------------------------------

# Read in and subset
upsig_knownmotifs <- read_delim("data/homer_upsig_deGenes_pval01_l2fc2/knownResults.txt") %>%
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
  mutate(motifLogo = paste0("data/homer_upsig_deGenes_pval01_l2fc2/knownResults/known", row_number(), ".logo.svg"))
  

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

# Downregulated motifs ----------------------------------------------------

# Read in and subset for top 3
downsig_knownmotifs <- read_delim("data/homer_downsig_deGenes_pval01_l2fc2/knownResults.txt") %>%
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
  mutate(motifLogo = paste0("data/homer_downsig_deGenes_pval01_l2fc2/knownResults/known", row_number(), ".logo.svg"))

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

grDevices::cairo_pdf(paste0("plots/tfmotif_expression.pdf"), 
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

