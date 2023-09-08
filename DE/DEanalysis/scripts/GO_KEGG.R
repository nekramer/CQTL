library(tidyverse)
library(rrvgo)
library(RColorBrewer)
library(ggnewscale)
library(patchwork)
source("scripts/plotting_utils.R")

# GO ----------------------------------------------------------------------

# Get reduced, significant GO terms for each category
upsig_go_data <- 
  read_delim("data/homer_upsig_deGenes_pval01_l2fc2/biological_process.txt") %>%
  mutate(pval = exp(1)^logP) %>%
  filter(pval < 0.01)
upsig_go <- reduceGO(upsig_go_data,
         category = "Upregulated")


downsig_go_data <- 
  read_delim("data/homer_downsig_deGenes_pval01_l2fc2/biological_process.txt") %>%
  mutate(pval = exp(1)^logP) %>%
  filter(pval < 0.01)
downsig_go <- reduceGO(downsig_go_data,
         category = "Downregulated")


# Select 5 each for plotting
upsig_go_plotting <- upsig_go %>%
  filter(parentTerm %in% c("response to cytokine", "cell surface receptor signaling pathway",
                           "collagen catabolic process", "regulation of cell-cell adhesion",
                           "acute inflammatory response")) %>%
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

ggsave(filename = "plots/GO_barplots.pdf", width = 5, height = 8, units = "in")


# KEGG --------------------------------------------------------------------

# Plot top 20 significant for each category
upsig_kegg_data <- read_delim("data/homer_upsig_deGenes_pval01_l2fc2/kegg.txt") %>%
  mutate(pval = exp(1)^logP) %>%
  filter(pval < 0.01) %>%
  distinct(Term, .keep_all = TRUE) %>%
  mutate(`-log10pval` = -log10(pval)) %>%
  mutate(category = "Upregulated")
downsig_kegg_data <- read_delim("data/homer_downsig_deGenes_pval01_l2fc2/kegg.txt") %>%
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

ggsave(filename = "plots/KEGG_barplots.pdf", width = 5, height = 8, units = "in")


# TF MOTIFS AND TF GENE EXPRESSION  ---------------------------------------




