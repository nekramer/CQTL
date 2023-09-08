library(tidyverse)
library(DESeq2)
library(ggrepel)
source("scripts/plotting_utils.R")

load('../processing/output/2023-04-19_gse.rda')
load("data/differential_expression_dds.rda")

# MA plot -----------------------------------------------------------------

# Prepare data - lfc shrink and grab data from plotMA
de_genes_shrink <- lfcShrink(dds, coef = "Condition_FNF_vs_CTL")
MAplot_data <- plotMA(de_genes_shrink, alpha = 0.01, returnData = TRUE) %>%
  rownames_to_column(var = "gene_id") %>%
  # Add column for up/down/no DE category
  mutate(category = ifelse(isDE == TRUE, 
                           ifelse(lfc > 0, "upDE", "downDE"), 
                           "noDE")) %>%
  # Get gene symbols for labeling
  left_join(as.data.frame(rowData(gse)) %>% 
              dplyr::select(c("gene_id", "symbol")), by = "gene_id") %>%
  # Labels for genes with strongest lfc
  mutate(symbol = ifelse(abs(lfc) > 5, symbol, ""))


ma_plot <- ggplot(MAplot_data, mapping = aes(x = mean, y = lfc, color = category)) + 
  geom_point(size = 0.75) +
  geom_text_repel(aes(label = symbol, family = "Helvetica",
                                           segment.color = "transparent"),
                  position = position_nudge_repel(x = 0.1, y = 0.01),
                  box.padding = 0.5,
                  size = 3,
                  max.overlaps = 5) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey") +
  geom_hline(yintercept = -1, linetype = "dashed", color = "grey") +
  scale_x_continuous(trans = "log10", 
                     breaks = c(1, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7)) +
  scale_color_manual(values = c("#78A1Cd", "grey", "#FBBE67")) +
  scale_y_continuous(breaks = seq(-20, 20, 5)) +
  theme_custom_scatterplot() +
  theme(legend.position = "none") +
  labs(x = "Mean of Normalized Counts",
       y = "Log Fold Change") 

ggsave(filename = "plots/ma_plot.pdf", width = 10, height = 11, units = "in")
save(ma_plot, file = "data/ma_plot.rda")


# Volcano plot ------------------------------------------------------------

de_genes_volcano <- as.data.frame(de_genes_shrink) %>%
  rownames_to_column(var = "gene_id") %>%
  # Set any padj that were 0 to smallest number to avoid Inf
  mutate(padj = ifelse(padj == 0, 5e-324, padj)) %>%
  # Calculate -log10pval for plotting
  mutate(`-log10pval` = -log10(padj)) %>%
  # Add DE type category for coloring
  mutate(category = ifelse(padj < 0.01, 
                           ifelse(log2FoldChange > 0, "upDE", "downDE"), 
                           "noDE")) %>%
  # Get gene symbols for labeling
  left_join(as.data.frame(rowData(gse)) %>% 
              dplyr::select(c("gene_id", "symbol")), by = "gene_id") %>%
  distinct(symbol, .keep_all = TRUE) %>%
  # Labels for genes with strongest lfc
  mutate(symbol = ifelse(abs(log2FoldChange) > 5, symbol, ""))


# Labels for ones with extremely large p-values
pval_labels <- de_genes_volcano %>% 
  filter(padj == 5e-324) %>%
  dplyr::select(-symbol) %>%
  left_join(as.data.frame(rowData(gse)) %>% 
              dplyr::select(c("gene_id", "symbol")), by = "gene_id")


  
volcano_plot <- ggplot(de_genes_volcano, aes(x = log2FoldChange, 
                                             y = `-log10pval`, 
                                             color = category)) +
  geom_point(size = 0.75) +
  geom_text_repel(aes(label = symbol, family = "Helvetica",
                      segment.color = "transparent"),
                  position = position_nudge_repel(y = 1),
                  box.padding = 0.5,
                  size = 3,
                  max.overlaps = 5) +
  geom_text_repel(data = pval_labels, aes(label = symbol, family = "Helvetica"),
                  position = position_nudge_repel(y = 15),
                  hjust = 0.5,
                  direction = "x",
                  force = 0.275,
                  box.padding = 0.5,
                  size = 3,
                  max.overlaps = 50,
                  segment.size = unit(0.25, "mm"),
                  segment.angle = 45) +
  scale_x_continuous(breaks = seq(-20, 20, 5), limits = c(-25, 20)) +
  scale_y_continuous(breaks = seq(0, 350, 50), limits = c(0, 350)) + 
  scale_color_manual(values = c("#78A1Cd", "grey", "#FBBE67")) +
  theme_custom_scatterplot() +
  theme(legend.position = "none") +
  labs(x = "Log Fold Change",
       y = "-log10(p-value)") 

ggsave(filename = "plots/volcano_plot.pdf", width = 16, height = 10, units = "in")
save(volcano_plot, file = "data/volcano_plot.rda")
