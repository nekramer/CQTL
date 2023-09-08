library(tidyverse)
library(patchwork)
library(ggtext)
library(colorspace)
library(grid)
library(VennDiagram)
library(ggpubr)
library(gghighlight)
library(ggrepel)
library(org.Hs.eg.db)
source("scripts/plotting_utils.R")

#### UNION OF SIG UNTREATED AND TREATED
# Union of sig untreated and treated genes --------------------------------
ctl_sig_genes <- read_csv("data/ctl_sexDE_pval01.csv")
fnf_sig_genes <- read_csv("data/fnf_sexDE_pval01.csv")

union_sig_genes <- union(ctl_sig_genes$gene_id, fnf_sig_genes$gene_id)

sex_genes <- bind_rows(ctl_sig_genes %>%
            filter(gene_id %in% union_sig_genes),
          fnf_sig_genes %>%
            filter(gene_id %in% union_sig_genes)) %>%
  distinct(gene_id, .keep_all = TRUE) %>%
  mutate(log2FC_dir = ifelse(log2FoldChange < 0, "female", "male"))

de_genes_results <- read_csv("data/de_genes_results.csv") %>%
  group_by(seqnames) %>%
  summarize(chrom_count = n())

# Chrom scatterplot with split log2fc  -----------------------------------------

# Get counts per chromosome and calculate avg abs log2fc of chrom
chrom_counts_l2fc_split <- sex_genes %>%
  mutate(`-log10pval` = -log10(padj)) %>%
  group_by(seqnames, log2FC_dir) %>%
  summarise(count = n(),
            avg_log2FoldChange = mean(log2FoldChange),
            med_log10pval = median(`-log10pval`)) %>%
  left_join(de_genes_results, by = "seqnames") %>%
  mutate(scaled_count = (count/chrom_count)*100) %>%
  mutate(alpha = ifelse(seqnames %in% c("X", "Y"), "yes", "no"))


chrom_counts_l2fc_split$chrom <- factor(chrom_counts_l2fc_split$seqnames,
                                        levels = c(as.character(1:22), "X", "Y"))


ggplot(chrom_counts_l2fc_split, aes(x = count, y = avg_log2FoldChange, 
                                    fill = log2FC_dir)) +
  geom_point(aes(size = med_log10pval, alpha = alpha), pch = 21, stroke = NA) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_text_repel(aes(label = ifelse(chrom %in% c("X", "Y"), 
                                     as.character(chrom), NA)),
                  nudge_x = -1, fontface = "bold", size = 5, 
                  min.segment.length = Inf,
                  color = "black") +
  scale_alpha_manual(values = c(0.7, 1)) +
  scale_fill_manual(values = c(sexColors[["F"]], sexColors[["M"]])) +
  scale_x_continuous(name = "Number of significant sex-specific genes",
                     limits = c(0, 22),
                     expand = c(0, 0)) +
  scale_y_continuous(name = "Average log~2~(fold change)",
                     limits = c(-2,8.2),
                     expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill = NA),
        panel.background = element_blank(),
        axis.line = element_line(),
        axis.ticks.x = element_blank(),
        axis.ticks.length.y = unit(-0.1, "cm"),
        axis.ticks.y = element_line(color = "black"),
        text = element_text(family = "Helvetica"),
        axis.text = element_text(color = "black", size = 11),
        axis.title.y = element_markdown(size = 12),
        axis.title.x = element_text(size = 12),
        legend.position = "none")

ggsave(filename = "plots/sex_chrom_splitl2fc_scatter.pdf", width = 6,
       height = 6, units = "in")


# x-axis with scaled percentage of sex-specific genes
ggplot(chrom_counts_l2fc_split, aes(x = scaled_count, y = avg_log2FoldChange, 
                                    fill = log2FC_dir)) +
  geom_point(aes(size = med_log10pval, alpha = alpha), pch = 21, stroke = NA) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_text_repel(aes(label = ifelse(chrom %in% c("X", "Y"), 
                                     as.character(chrom), NA)),
                  nudge_x = 1.5, fontface = "bold", size = 5, 
                  min.segment.length = Inf,
                  color = "black") +
  scale_alpha_manual(values = c(0.5, 1)) +
  scale_fill_manual(values = c(sexColors[["F"]], sexColors[["M"]])) +
  scale_x_continuous(name = "Percentage of sex-specific chromosome genes",
                     limits = c(0, 55),
                     expand = c(0, 0)) +
  scale_y_continuous(name = "Average log~2~(fold change)",
                     limits = c(-2,8.2),
                     expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill = NA),
        panel.background = element_blank(),
        axis.line = element_line(),
        axis.ticks.x = element_blank(),
        axis.ticks.length.y = unit(-0.1, "cm"),
        axis.ticks.y = element_line(color = "black"),
        text = element_text(family = "Helvetica"),
        axis.text = element_text(color = "black", size = 11),
        axis.title.y = element_markdown(size = 12),
        axis.title.x = element_text(size = 12),
        legend.position = "none")


ggsave(filename = "plots/sex_chrom_splitl2fc_scatter_scaled.pdf", width = 6,
       height = 6, units = "in")



# Chromosome counts -------------------------------------------------------


# Get counts per chromosome, per sex
chrom_counts <- sex_genes %>%
  group_by(seqnames, log2FC_dir) %>%
  summarise(count = n())
chrom_counts$chrom <- factor(chrom_counts$seqnames,
                             levels = c(as.character(1:22), "X", "Y"))

# Lollipop plot of counts
ggplot(chrom_counts, aes(x = chrom, y = count, fill = log2FC_dir, color = log2FC_dir)) +
  geom_linerange(aes(x = chrom, ymin = 0, ymax = count),
                 position = position_dodge(width = 0.5),
                 linewidth = 1) +
  geom_point(position = position_dodge(width = 0.5),
             size = 3) +
  scale_fill_manual(values = c(sexColors[["F"]], sexColors[["M"]])) +
  scale_color_manual(values = c(sexColors[["F"]], sexColors[["M"]])) +
  scale_x_discrete(name = "Chromosome", drop = FALSE) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 30),
                     labels = c(10, 20, 30),
                     breaks = c(10, 20, 30),
                     name = 'Number of 
                     <span style = "color:#4788BA;">**male**</span> and
                     <span style = "color:#DD8492;">**female**</span><br> significant genes') +
  theme(text = element_text(family = "Helvetica"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = "none",
        axis.ticks.x = element_blank(),
        axis.ticks.length.y = unit(-0.15, "cm"),
        axis.line = element_line(),
        axis.text = element_text(color = "black", size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_markdown(size = 12))

ggsave(filename = "plots/sex_de_pval01_union_sex_chromlollipops.pdf", 
       width = 9, height = 3.5, units = "in")


# Overlap with sex-specific genes from other tissues ----------------------

# Pull union list
chond_sbgenes <- sex_genes %>%
  mutate(sex = ifelse(log2FoldChange < 0, "female", "male")) %>%
  mutate(tissue = "Chondrocytes") %>%
  dplyr::rename(effsize = log2FoldChange) %>%
  dplyr::rename(effsize_se = lfcSE)



# chond_sbgenes <- read_csv("data/ctl_sex_shrink.csv") %>%
#   filter(gene_id %in% union_sig_genes) %>%
#   mutate(sex = ifelse(log2FoldChange < 0, "female", "male")) %>%
#   mutate(tissue = "Chondrocytes") %>%
#   dplyr::rename(effsize = log2FoldChange) %>%
#   dplyr::rename(effsize_se = lfcSE)

# Read in gtex sex-biased genes
gtex_signif_sbgenes <- 
  read_delim("data/GTEx_Analysis_v8_sbgenes/signif.sbgenes.txt") %>%
  # convert IDs to ones compatible with ours
  mutate(gene_id = gsub("\\..*", "", gene)) %>%
  # Flip sign of effsize to match ours (gtex positive is female and negative is male)
  mutate(effsize = -1*effsize) %>%
  mutate(sex = ifelse(effsize < 0, "female", "male")) %>%
  # filter for genes in chond_sbgenes
  filter(gene_id %in% chond_sbgenes$gene_id) %>%
  mutate(tissue = gsub("_", " ", tissue))
  
# # Get gene symbols
# gtex_geneSymbols <- AnnotationDbi::select(org.Hs.eg.db,
#                                           keys = gtex_signif_sbgenes$gene_id,
#                                           columns = c("ENSEMBL", "SYMBOL"),
#                                           keytype = "ENSEMBL")

# Join with our data
all_data_sbgenes <- bind_rows(chond_sbgenes %>% dplyr::select(gene_id, tissue, sex, effsize, effsize_se),
          gtex_signif_sbgenes %>% dplyr::select(gene_id, tissue, sex, effsize, effsize_se)) %>%
  complete(gene_id, tissue) %>% 
  left_join(chond_sbgenes %>% dplyr::select(gene_id, symbol), by = "gene_id")
  
# Order genes by how many overlap
geneOverlaps <- all_data_sbgenes %>%
  group_by(symbol) %>%
  filter(!is.na(effsize)) %>%
  summarize(nOverlap = n()) %>%
  arrange(nOverlap)

all_data_sbgenes <- left_join(all_data_sbgenes, geneOverlaps, by = "symbol") %>%
  group_by(nOverlap) %>%
  arrange(sex, .by_group = TRUE) %>%
  ungroup()


all_data_sbgenes$tissue <- factor(all_data_sbgenes$tissue,
                                  levels = c("Chondrocytes", unique(gtex_signif_sbgenes$tissue)))

all_data_sbgenes$symbol <- factor(all_data_sbgenes$symbol,
                                   levels = unique(all_data_sbgenes$symbol))

all_data_sbgenes <- all_data_sbgenes %>%
  mutate(xmin = 0.5, xmax = 45.5, 
         y_position = as.numeric(symbol),
         ymin = y_position - 0.5,
         ymax = y_position + 0.5) %>%
  pivot_longer(cols = c(xmin, xmax), values_to = "x", names_to="xmin_xmax") %>%
  dplyr::select(-xmin_xmax) %>%
  mutate(fill = ifelse(y_position %% 2 == 0, "a", "b"))

fontFaces <- all_data_sbgenes %>%
  reframe(face = ifelse(nOverlap == 1, "bold", "plain"), .by = symbol) %>%
  distinct()

# Heatmap of M vs F for each gene
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

# ggsave(filename = "plots/sex_de_celltypeOverlap.pdf",
#        width = 11, height = 14, units = "in")



chond_sbgenes <- left_join(chond_sbgenes, geneOverlaps, by = "symbol") %>%
  arrange(nOverlap)
chond_sbgenes$symbol <- factor(chond_sbgenes$symbol,
                                  levels = unique(all_data_sbgenes$symbol))


chond_sbgenes <- chond_sbgenes %>%
  mutate(xmin = -5, xmax = 20, 
         y_position = as.numeric(symbol),
         ymin = y_position - 0.5,
         ymax = y_position + 0.5) %>%
  pivot_longer(cols = c(xmin, xmax), values_to = "x", names_to="xmin_xmax") %>%
  dplyr::select(-xmin_xmax) %>%
  mutate(fill = ifelse(y_position %% 2 == 0, "a", "b"))



celloverlap_chondl2fc <- ggplot(chond_sbgenes, aes(x = effsize, y = y_position)) +
  geom_ribbon(aes(x = x, ymin = ymin, ymax = ymax, group = y_position, fill = fill), 
              inherit.aes = FALSE, alpha = 0.2) +
  scale_fill_manual(values = c("grey50", "grey75")) +
  #ggnewscale::new_scale_fill() +
  geom_vline(xintercept = 0, lty = 2) +
  geom_segment(aes(x = effsize-1.96*effsize_se, xend = effsize+1.96*effsize_se,
                   y = y_position, yend = y_position), color = "grey25") +
  geom_point(aes(color = sex)) +
  coord_cartesian(clip = "off") +
  scale_color_manual(values = c(sexColors[["F"]], sexColors[["M"]])) +
  scale_x_continuous(limits = c(-5, 20), name = "Gene sex log~2~(fold change)<br> in chondrocytes",
                     expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0))+
  theme(axis.title.x = element_markdown(family = "Helvetica"),
        axis.title.y = element_blank(),
        legend.position = "none",
        panel.border = element_rect(fill = NA, color = "grey25"),
        panel.background = element_blank(),
        axis.text.x = element_text(color = "black", size = 10),
        axis.text.y = element_blank(),
        axis.line.x = element_line(color = "grey25"),
        axis.ticks = element_blank(),
        text = element_text(family = "Helvetica"), 
        plot.margin = margin(0, 20, 0, 0))

# ggsave(filename = "plots/sex_de_celltypeOverlap_chondl2fc.pdf",
#        width = 3, height = 14, units = "in")

# 
sex_cell_heatmap + celloverlap_chondl2fc +
  plot_layout(ncol = 2, widths = c(4, 1))


ggsave(filename = "plots/sex_celltypeOverlap_heatmap_l2fc.pdf",
       width = 13, height = 15, units = "in")


#### PBS AND FNF SEPARATED

# PBS and FNF separated ---------------------------------------------------
ctl_sig_genes <- read_csv("data/ctl_sexDE_pval01.csv") %>%
  mutate(condition = "PBS") %>%
  mutate(sex = ifelse(log2FoldChange < 0, "Female", "Male"))
fnf_sig_genes <- read_csv("data/fnf_sexDE_pval01.csv") %>%
  mutate(condition = "FN-f") %>%
  mutate(sex = ifelse(log2FoldChange < 0, "Female", "Male"))

all_sig_sex_genes <- bind_rows(ctl_sig_genes,
                               fnf_sig_genes)
all_sig_sex_genes$seqnames <- factor(all_sig_sex_genes$seqnames,
                                     levels = c(as.character(1:22), "X", "Y"))
# Stacked bar plot of chrom counts ----------------------------------------

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


ggsave(filename = "plots/sex_fnfpbs_numGenes_barplot.pdf",
       width = 4, height = 3.5, units = "in")


# Lollipop plot of chrom counts -------------------------------------------

ctl_desex_pval01 <- read_csv("data/ctl_sexDE_pval01.csv") %>%
  mutate(padj = ifelse(padj == 0, .Machine$double.xmin, padj)) %>%
  mutate(condition = "Control") %>%
  mutate(id = "pval")

fnf_desex_pval01 <- read_csv("data/fnf_sexDE_pval01.csv") %>%
  mutate(padj = ifelse(padj == 0, .Machine$double.xmin, padj)) %>%
  mutate(condition = "FN-f") %>%
  mutate(id = "pval")


ctl_fnf_desex_pval01 <- bind_rows(ctl_desex_pval01,
                                  fnf_desex_pval01)
ctl_fnf_desex_pval01$chrom <- factor(ctl_fnf_desex_pval01$seqnames,
                                     levels = c(as.character(1:22), "X", "Y"))
ctl_fnf_desex_pval01 <- ctl_fnf_desex_pval01 %>%
  mutate(aut_sex = ifelse(chrom %in% c("X", "Y"), "sex", "aut")) %>%
  mutate(aut_sex2 = ifelse(chrom %in% c("X", "Y"), seqnames, "aut"))



## Separate all chromosomes

ctl_fnf_desex_pval01_chrom_counts <- ctl_fnf_desex_pval01 %>%
  group_by(chrom, condition) %>%
  summarise(count = n())

ggplot(ctl_fnf_desex_pval01_chrom_counts, aes(x = chrom, y = count, fill = condition, color = condition)) +
  geom_linerange(aes(x = chrom, ymin = 0, ymax = count),
                 position = position_dodge(width = 0.5)) +
  geom_point(position = position_dodge(width = 0.5)) +
  scale_fill_manual(values = c("#EDB361", "#6996C7")) +
  scale_color_manual(values = c("#EDB361", "#6996C7")) +
  scale_x_discrete(name = "Chromosome", drop = FALSE) +
  scale_y_continuous(expand = c(0,0), name = "Count", limits = c(0, 30)) +
  theme_custom_general() +
  theme(text = element_text(family = "Helvetica"),
        legend.position = "none",
        plot.title = ggtext::element_markdown(size = 12)) +
  labs(title = 'Number of sex-specific differential genes in
       <span style = "color:#EDB361;">**control**</span> and
       <span style = "color:#6996C7;">**FN-f**</span>')

ggsave(filename = "plots/sex_de_pval01_chromlollipops.pdf", 
       width = 11, height = 6, units = "in")


## Autosomes vs sex chromosomes 
ctl_fnf_desex_pval01_autsex_counts <- ctl_fnf_desex_pval01 %>%
  group_by(aut_sex, condition) %>%
  summarise(count = n())

ggplot(ctl_fnf_desex_pval01_autsex_counts, aes(x = aut_sex, y = count, fill = condition, color = condition)) +
  geom_linerange(aes(x = aut_sex, ymin = 0, ymax = count),
                 position = position_dodge(width = 0.5)) +
  geom_point(position = position_dodge(width = 0.5)) +
  scale_fill_manual(values = c("#EDB361", "#6996C7")) +
  scale_color_manual(values = c("#EDB361", "#6996C7")) +
  scale_x_discrete(labels = c("Autosomes", "Sex Chromosomes")) +
  scale_y_continuous(expand = c(0,0), name = "Count", limits = c(0, 55)) +
  theme_custom_general() +
  theme(text = element_text(family = "Helvetica"),
        legend.position = "none",
        axis.title.x = element_blank(), 
        plot.title = ggtext::element_markdown(size = 12)) +
  labs(title = 'Number of sex-specific differential genes in
       <span style = "color:#EDB361;">**control**</span> and
       <span style = "color:#6996C7;">**FN-f**</span>')

ggsave(filename = "plots/sex_de_pval01_autosome_sexchromlollipops.pdf", 
       width = 6, height = 6, units = "in")



## Autosomes vs sex chromosomes  v2
ctl_fnf_desex_pval01_autsex2_counts <- ctl_fnf_desex_pval01 %>%
  group_by(aut_sex2, condition) %>%
  summarise(count = n())

ggplot(ctl_fnf_desex_pval01_autsex2_counts, aes(x = aut_sex2, y = count, 
                                                fill = condition, 
                                                color = condition)) +
  geom_linerange(aes(x = aut_sex2, ymin = 0, ymax = count),
                 position = position_dodge(width = 0.5)) +
  geom_point(position = position_dodge(width = 0.5)) +
  scale_fill_manual(values = c("#EDB361", "#6996C7")) +
  scale_color_manual(values = c("#EDB361", "#6996C7")) +
  scale_x_discrete(labels = c("Autosomes", "X Chromosome", "Y Chromosome")) +
  scale_y_continuous(expand = c(0,0), name = "Count", limits = c(0, 55)) +
  theme_custom_general() +
  theme(text = element_text(family = "Helvetica"),
        legend.position = "none",
        axis.title.x = element_blank(), 
        plot.title = ggtext::element_markdown(size = 12)) +
  labs(title = 'Number of sex-specific differential genes in
       <span style = "color:#EDB361;">**control**</span> and
       <span style = "color:#6996C7;">**FN-f**</span>')

ggsave(filename = "plots/sex_de_pval01_autosome_sexchromlollipops_v2.pdf", 
       width = 6, height = 6, units = "in")

 # Chrom p-value boxplots -------------------------------------------------------

ctl_desex_pval01 <- read_csv("data/ctl_sexDE_pval01.csv") %>%
  mutate(padj = ifelse(padj == 0, .Machine$double.xmin, padj)) %>%
  mutate(condition = "Control") %>%
  mutate(id = "pval") %>%
  mutate(log2FC_dir = ifelse(log2FoldChange < 0, "-", "+"))

fnf_desex_pval01 <- read_csv("data/fnf_sexDE_pval01.csv") %>%
  mutate(padj = ifelse(padj == 0, .Machine$double.xmin, padj)) %>%
  mutate(condition = "FN-f") %>%
  mutate(id = "pval") %>%
  mutate(log2FC_dir = ifelse(log2FoldChange < 0, "-", "+"))


ctl_fnf_desex_pval01 <- bind_rows(ctl_desex_pval01,
                                  fnf_desex_pval01)
ctl_fnf_desex_pval01$chrom <- factor(ctl_fnf_desex_pval01$seqnames,
                                     levels = c(as.character(1:22), "X", "Y"))

ctl_fnf_desex_pval01 <- ctl_fnf_desex_pval01 %>%
  mutate(aut_sex = ifelse(chrom %in% c("X", "Y"), "sex", "aut")) %>%
  mutate(aut_sex2 = ifelse(chrom %in% c("X", "Y"), seqnames, "aut")) %>%
  mutate(log10pval = -log10(padj))


## Separate all chromosomes
ggplot(ctl_fnf_desex_pval01, aes(x = chrom, y = -log10(padj), fill = condition, color = condition)) +
  geom_boxplot(position = position_dodge(width = 1)) +
  scale_fill_manual(values = c("#EDB361", "#6996C7")) +
  scale_color_manual(values = c("#EDB361", "#6996C7")) +
  scale_x_discrete(name = "Chromosome", drop = FALSE) +
  scale_y_continuous(expand = c(0,0), 
                     name = "-log~10~(padj)") +
  theme_custom_general() +
  theme(text = element_text(family = "Helvetica"),
        legend.position = "none",
        plot.title = ggtext::element_markdown(size = 12),
        axis.title.y = element_markdown()) +
  labs(title = 'Significance of sex-specific differential genes in
       <span style = "color:#EDB361;">**control**</span> and
       <span style = "color:#6996C7;">**FN-f**</span>')

ggsave(filename = "plots/sex_de_pval01_chromsignificance.pdf", 
       width = 11, height = 6, units = "in")


## Autosomes vs sex chromosomes 
ggplot(ctl_fnf_desex_pval01, aes(x = aut_sex, y = -log10(padj), fill = condition, color = condition)) +
  geom_boxplot(position = position_dodge(width = 1)) +
  scale_fill_manual(values = c("#EDB361", "#6996C7")) +
  scale_color_manual(values = c("#EDB361", "#6996C7")) +
  scale_x_discrete(labels = c("Autosomes", "Sex Chromosomes")) +
  scale_y_continuous(expand = c(0,0), 
                     name = "-log~10~(padj)",
                     limits = c(0, 325)) +
  theme_custom_general() +
  theme(text = element_text(family = "Helvetica"),
        legend.position = "none",
        axis.title.x = element_blank(),
        plot.title = ggtext::element_markdown(size = 12),
        axis.title.y = element_markdown()) +
  labs(title = 'Significance of sex-specific differential genes in
       <span style = "color:#EDB361;">**control**</span> and
       <span style = "color:#6996C7;">**FN-f**</span>')

ggsave(filename = "plots/sex_de_pval01_autosome_sexchromsignificance.pdf", 
       width = 6, height = 6, units = "in")


## Autosomes vs sex chromosomes v2
ggplot(ctl_fnf_desex_pval01 %>%
         dplyr::select(gene_id, aut_sex2, log10pval, condition), 
       aes(x = aut_sex2, y = log10pval,
           fill = condition, color = condition)) +
  geom_boxplot(position = position_dodge(width = 1)) +
  stat_compare_means(aes(group = condition), comparisons = list(c("aut", "X"), c("aut", "Y")),
                     label = "p.signif",
                     tip.length = 0,
                     show.legend = TRUE) +
  scale_fill_manual(values = c("#EDB361", "#6996C7")) +
  scale_color_manual(values = c("#EDB361", "#6996C7")) +
  scale_x_discrete(labels = c("Autosomes", "X Chromosome", "Y Chromosome")) +
  scale_y_continuous(expand = c(0,0), 
                     name = "-log~10~(padj)",
                     limits = c(0, 400)) +
  theme_custom_general() +
  theme(text = element_text(family = "Helvetica"),
        legend.position = "none",
        axis.title.x = element_blank(),
        plot.title = ggtext::element_markdown(size = 12),
        axis.title.y = element_markdown()) +
  labs(title = 'Significance of sex-specific differential genes in
       <span style = "color:#EDB361;">**control**</span> and
       <span style = "color:#6996C7;">**FN-f**</span>')
# ggsave(filename = "plots/sex_de_pval01_autosome_sexchromsignificance_v2.pdf", 
#        width = 6, height = 6, units = "in")
#****: p <= 0.0001
ggsave(filename = "plots/sex_de_pval01_autosome_sexchromsignificance_boxsignif.pdf",
       width = 6, height = 6, units = "in")

## Autosomes vs sex chromosomes v3
ggplot(ctl_fnf_desex_pval01, 
       aes(x = aut_sex2, y = log10pval,
           fill = condition, color = condition)) +
  geom_boxplot(position = position_dodge(width = 1)) +
  scale_fill_manual(values = c("#EDB361", "#6996C7")) +
  scale_color_manual(values = c("#EDB361", "#6996C7")) +
  facet_wrap(vars(log2FC_dir), nrow = 2, labeller = as_labeller(c("-" = "Female", "+" = "Male"))) +
  scale_x_discrete(labels = c("Autosomes", "X Chromosome", "Y Chromosome")) +
  scale_y_continuous(expand = c(0,0), 
                     name = "-log~10~(padj)",
                     limits = c(0, 400)) +
  theme_custom_general() +
  theme(text = element_text(family = "Helvetica"),
        legend.position = "none",
        axis.title.x = element_blank(),
        plot.title = ggtext::element_markdown(size = 12),
        axis.title.y = element_markdown(),
        strip.background = element_blank(),
        strip.text = element_text(family = "Helvetica", size = 12)) +
  labs(title = 'Significance of sex-specific differential genes in
       <span style = "color:#EDB361;">**control**</span> and
       <span style = "color:#6996C7;">**FN-f**</span>')


ggsave(filename = "plots/sex_de_pval01_autosome_sexchromsignificance_malefemalesplit.pdf", 
       width = 6, height = 10, units = "in")
# All p-values (not just sig) --------------------------------------------------
load("data/dds_sex_ctl.rda")
load("data/dds_sex_fnf.rda")

ctl_sexgenes_shrink <- lfcShrink(dds_sex_ctl, coef = "Sex_M_vs_F",
                             format = "GRanges") %>%
  plyranges::names_to_column("gene_id") %>%
  as_tibble() %>%
  mutate(condition = "Control")

fnf_sexgenes_shrink <- lfcShrink(dds_sex_fnf, coef = "Sex_M_vs_F",
                                 format = "GRanges") %>%
  plyranges::names_to_column("gene_id") %>%
  as_tibble() %>%
  mutate(condition = "FN-f")

all_sexgenes_shrink <- bind_rows(ctl_sexgenes_shrink, fnf_sexgenes_shrink)
all_sexgenes_shrink$chrom <- factor(all_sexgenes_shrink$seqnames,
                                     levels = c(as.character(1:22), "X", "Y"))
all_sexgenes_shrink <- all_sexgenes_shrink %>%
  mutate(padj = ifelse(padj == 0, .Machine$double.xmin, padj)) %>%
  mutate(aut_sex2 = ifelse(chrom %in% c("X", "Y"), as.character(chrom), "aut")) %>%
  mutate(log10pval = -log10(padj)) %>%
  mutate(sig = ifelse(padj < 0.01, "sig", "notsig")) %>%
  na.omit() %>%
  mutate(cond_autsex2 = paste0(condition, "_", aut_sex2)) %>%
  mutate(cond_autsex2_sig = paste0(condition, "_", aut_sex2, "_", sig))

stat_means <- compare_means(log10pval ~ sig,
                            all_sexgenes_shrink,
                            group.by = c("condition", "aut_sex2"))




ggplot(all_sexgenes_shrink, aes(x = aut_sex2, y = log10pval, fill = sig, color = sig)) +
  facet_wrap(vars(condition)) +
  geom_boxplot(position = position_dodge(width = 1)) +
  scale_fill_manual(values = c("grey75", "dodgerblue")) +
  scale_color_manual(values = c("grey75", "dodgerblue")) +
  stat_compare_means(label = "p.signif") +
  scale_x_discrete(labels = c("Autosomes", "X Chromosome", "Y Chromosome")) +
  scale_y_continuous(name = "-log~10~(padj)",
                     limits = c(0, 350)) +
  theme_custom_general() +
  theme(text = element_text(family = "Helvetica"),
        legend.position = "none",
        axis.title.x = element_blank(),
        plot.title = ggtext::element_markdown(size = 12),
        axis.title.y = element_markdown(),
        strip.background = element_blank()) +
  labs(title = 'Comparison of <span style = "color:grey75;">non-significant</span> and 
      <span style = "color:dodgerblue;">significant</span> sex-specific genes in
       control and FN-f')

ggsave(filename = "plots/sex_nonsig_vs_sig_pvals.pdf",
       width = 8, height = 6, units = "in")

# LFC distributions by chromosome -----------------------------------------------

ctl_desex_pval01 <- read_csv("data/ctl_sexDE_pval01.csv") %>%
  mutate(padj = ifelse(padj == 0, .Machine$double.xmin, padj)) %>%
  mutate(condition = "Control") %>%
  mutate(log2FC_dir = ifelse(log2FoldChange < 0, "-", "+"))

fnf_desex_pval01 <- read_csv("data/fnf_sexDE_pval01.csv") %>%
  mutate(padj = ifelse(padj == 0, .Machine$double.xmin, padj)) %>%
  mutate(condition = "FN-f") %>%
  mutate(log2FC_dir = ifelse(log2FoldChange < 0, "-", "+"))
ctl_fnf_desex_pval01 <- bind_rows(ctl_desex_pval01,
                                  fnf_desex_pval01)
ctl_fnf_desex_pval01$chrom <- factor(ctl_fnf_desex_pval01$seqnames,
                                     levels = c(as.character(1:22), "X", "Y"))

ctl_fnf_desex_pval01 <- ctl_fnf_desex_pval01 %>%
  mutate(aut_sex = ifelse(chrom %in% c("X", "Y"), "sex", "aut")) %>%
  mutate(aut_sex2 = ifelse(chrom %in% c("X", "Y"), seqnames, "aut"))

## Separate all chromosomes
ggplot(ctl_fnf_desex_pval01, aes(x = chrom, y = abs(log2FoldChange), fill = condition, color = condition)) +
  geom_boxplot(position = position_dodge(width = 1)) +
  scale_fill_manual(values = c("#EDB361", "#6996C7")) +
  scale_color_manual(values =  c("#EDB361", "#6996C7")) +
  facet_wrap(vars(log2FC_dir), nrow = 2, labeller = as_labeller(c("-" = "Female", "+" = "Male"))) +
  scale_x_discrete(name = "Chromosome", drop = FALSE) +
  scale_y_continuous(expand = c(0,0), 
                     name = "log~2~ Fold Change",
                     limits = c(0, 12),
                     breaks = seq(0, 12, 2)) +
  theme_custom_general() +
  theme(text = element_text(family = "Helvetica"),
        legend.position = "none",
        plot.title = ggtext::element_markdown(size = 12),
        axis.title.y = element_markdown(),
        strip.background = element_blank(),
        strip.text = element_text(family = "Helvetica", size = 12)) +
  labs(title = 'Effect size of sex-specific differential genes in
       <span style = "color:#EDB361;">**control**</span> and
       <span style = "color:#6996C7;">**FN-f**</span>')

ggsave(filename = "plots/sex_de_pval01_chrom_effectsizes.pdf", 
       width = 11, height = 10, units = "in")

## Autosomes vs sex chromosomes 
ggplot(ctl_fnf_desex_pval01, aes(x = aut_sex, y = abs(log2FoldChange), fill = condition, color = condition)) +
  geom_boxplot(position = position_dodge(width = 1)) +
  scale_fill_manual(values = c("#EDB361", "#6996C7")) +
  scale_color_manual(values =  c("#EDB361", "#6996C7")) +
  facet_wrap(vars(log2FC_dir), nrow = 2, labeller = as_labeller(c("-" = "Female", "+" = "Male"))) +
  scale_x_discrete(labels = c("Autosomes", "Sex Chromosomes")) +
  scale_y_continuous(expand = c(0,0), 
                     name = "log~2~ Fold Change",
                     limits = c(0, 12),
                     breaks = seq(0, 12, 2)) +
  theme_custom_general() +
  theme(text = element_text(family = "Helvetica"),
        legend.position = "none",
        axis.title.x = element_blank(),
        plot.title = ggtext::element_markdown(size = 12),
        axis.title.y = element_markdown(),
        strip.background = element_blank(),
        strip.text = element_text(family = "Helvetica", size = 12)) +
  labs(title = 'Effect size of sex-specific differential genes in
       <span style = "color:#EDB361;">**control**</span> and
       <span style = "color:#6996C7;">**FN-f**</span>')

ggsave(filename = "plots/sex_de_pval01_autosome_sexchrom_effectsizes.pdf", 
       width = 6, height = 10, units = "in")


compare_means(abs("log2FoldChange") ~ aut_sex2, ctl_fnf_desex_pval01)

## Autosomes vs sex chromosomes  v2
ggplot(ctl_fnf_desex_pval01, aes(x = aut_sex2, y = abs(log2FoldChange), fill = condition, color = condition)) +
  geom_boxplot(position = position_dodge(width = 1)) +
  stat_compare_means(aes(group = condition), comparisons = list(c("aut", "X"), c("aut", "Y")),
                     label = "p.signif",
                     tip.length = 0,
                     show.legend = TRUE) +
  scale_fill_manual(values = c("#EDB361", "#6996C7")) +
  scale_color_manual(values =  c("#EDB361", "#6996C7")) +
  facet_wrap(vars(log2FC_dir), nrow = 2, labeller = as_labeller(c("-" = "Female", "+" = "Male"))) +
  scale_x_discrete(labels = c("Autosomes", "X Chromosome", "Y Chromosome")) +
  scale_y_continuous(expand = c(0,0), 
                     name = "log~2~ Fold Change") +
  theme_custom_general() +
  theme(text = element_text(family = "Helvetica"),
        legend.position = "none",
        axis.title.x = element_blank(),
        plot.title = ggtext::element_markdown(size = 12),
        axis.title.y = element_markdown(),
        strip.background = element_blank(),
        strip.text = element_text(family = "Helvetica", size = 12)) +
  labs(title = 'Effect size of sex-specific differential genes in
       <span style = "color:#EDB361;">**control**</span> and
       <span style = "color:#6996C7;">**FN-f**</span>')


ggsave(filename = "plots/sex_de_pval01_autosome_sexchrom_effectsizes_v2.pdf", 
       width = 6, height = 10, units = "in")

# LFC summaries by condition ----------------------------------------------

ctl_desex_pval01 <- read_csv("data/ctl_sexDE_pval01.csv") %>%
  mutate(padj = ifelse(padj == 0, .Machine$double.xmin, padj)) %>%
  mutate(condition = "Control") %>%
  mutate(log2FC_dir = ifelse(log2FoldChange < 0, "-", "+"))

fnf_desex_pval01 <- read_csv("data/fnf_sexDE_pval01.csv") %>%
  mutate(padj = ifelse(padj == 0, .Machine$double.xmin, padj)) %>%
  mutate(condition = "FN-f") %>%
  mutate(log2FC_dir = ifelse(log2FoldChange < 0, "-", "+"))
ctl_fnf_desex_pval01 <- bind_rows(ctl_desex_pval01,
                                  fnf_desex_pval01)
ctl_fnf_desex_pval01$chrom <- factor(ctl_fnf_desex_pval01$seqnames,
                                     levels = c(as.character(1:22), "X", "Y"))
ctl_fnf_desex_pval01$log2FC_dir <- factor(ctl_fnf_desex_pval01$log2FC_dir, levels = c("-", "+"))

ggplot(ctl_fnf_desex_pval01, aes(x = condition, y = abs(log2FoldChange), fill = log2FC_dir, color = log2FC_dir)) +
  geom_boxplot(position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = c(sexColors[[2]], sexColors[[1]])) +
  scale_color_manual(values =  c(sexColors[[2]], sexColors[[1]])) +
  scale_y_continuous(expand = c(0,0), 
                     name = "log~2~ Fold Change") +
  theme_custom_general() +
  theme(text = element_text(family = "Helvetica"),
        legend.position = "none",
        plot.title = ggtext::element_markdown(size = 12),
        axis.title.y = element_markdown(),
        strip.background = element_blank(),
        strip.text = element_text(family = "Helvetica", size = 12),
        axis.title.x = element_blank()) +
  labs(title = 'Effect sizes of <span style = "color:#6e82b7;">**male**</span> and
  <span style = "color:#f492a5;">**female**</span> differential genes')

ggsave(filename = "plots/sex_de_pval01_condition_effectsizes.pdf", 
       width = 6, height = 6, units = "in")


# Sex-specific response gene slope plots ----------------------------------

plotSexResponseGene <- function(geneRow, dds, suffix = ""){
  
  # Gene ID
  geneid <- geneRow[6]
  # Gene symbol
  genesymbol <- geneRow[13]
  # MvsF Contrast stats
  padj <- format(as.numeric(geneRow[12]), digits = 3, scientific = TRUE)
  lfc <- round(as.numeric(geneRow[8]), digits = 3)
  
  # F CTLvsFNF stats
  F_stats <- as.data.frame(results(dds, name = "SexF.ConditionFNF")) %>% 
    rownames_to_column(var = "gene_id") %>%
    filter(gene_id == geneid)
  F_padj <- format(F_stats %>% pull(padj), digits = 3, scientific = TRUE)
  F_lfc <- round(F_stats %>% pull(log2FoldChange), digits = 3)
  
  F_df <- data.frame("Sex" = "F",
                     "count" = 100)
  
  # M CTLvsFNF stats
  M_stats <- as.data.frame(results(dds, name = "SexM.ConditionFNF")) %>% 
    rownames_to_column(var = "gene_id") %>%
    filter(gene_id == geneid)
  M_padj <- format(M_stats %>% pull(padj), digits = 3, scientific = TRUE)
  M_lfc <- round(M_stats %>% pull(log2FoldChange), digits = 3)
  M_df <- data.frame("Sex" = "M",
                     "count" = 100)
  
  
  # Grab count data for gene from plotCounts
  geneData <- plotCounts(dds = dds, gene = geneid, intgroup = c("Sex", "Condition"),
                         returnData = TRUE) %>%
    mutate(group = paste0(Sex, "_", Condition))
  
  
  # Plot with ggplot
  plot <- ggplot(geneData, aes(x = Condition, y = count, color = Sex)) +
    geom_jitter(width = 0.1, size = 2) +
    facet_wrap(~Sex, strip.position = "bottom") +
    scale_y_continuous(trans = "log2", breaks = c(100, 200, 500, 1000, 2000), 
                       name = "Normalized counts") +
    scale_x_discrete(labels = c("CTL", "FNF", "CTL", "FNF")) +
    scale_color_manual(values = c("#C06E8B", "#6F8CC7")) +
    geom_richtext(data = F_df, aes(x = 1.5, family = "Helvetica"),
                  size = 3,
                  label = paste0("padj = ", F_padj, "<br>", "log~2~FC = ", F_lfc),
                  fill = NA, label.color = NA) +
    geom_richtext(data = M_df, aes(x = 1.5, family = "Helvetica"),
                  size = 3,
                  label = paste0("padj = ", M_padj, "<br>", "log~2~FC = ", M_lfc),
                  fill = NA, label.color = NA) +
    theme_custom_scatterplot() +
    theme(legend.position = "none",
          axis.title.x = element_blank(),
          strip.placement = "outside",
          strip.background = element_blank(),
          strip.text = element_text(size = 10),
          strip.switch.pad.wrap = unit(0, "in"),
          plot.subtitle = element_markdown(hjust = 0.05, size = 10, color = "grey25"),
          plot.title = element_text(hjust = 0.05)) +
    ggtitle(label = genesymbol, subtitle = paste0("padj = ", padj, "<br>",
                                                  "log~2~FC = ", lfc))
  
  # Save
  ggsave(plot, filename = paste0("plots/sexspecifictreatment_",
                                 genesymbol, suffix, ".pdf"),
         width = 6, height = 6, units = "in" )
  return(plot)
  
} 

load("data/dds_sex_treatmentresponse.rda")
treatmentresponse_sexDE_genes <- read_csv("data/sexDE_treatmenteffect_pval01.csv")
apply(treatmentresponse_sexDE_genes, 1, plotSexResponseGene, dds = dds_sex_treatmentresponse)


load("data/dds_sex_treatmentresponse_moregenefilter.rda")
treatmentresponse_sexDE_genes_moregenefilter <- 
  read_csv("data/sexDE_treatmenteffect_moregenefilter_pval01.csv")

apply(treatmentresponse_sexDE_genes_moregenefilter, 
      1, plotSexResponseGene, dds = dds_sex_treatmentresponse_moregenefilter,
      suffix = "_genefilter")




# Venn diagram between ctl and fnf  -------------------------------------------

ctl_desex_pval01 <- read_csv("data/ctl_sexDE_pval01.csv") 

fnf_desex_pval01 <- read_csv("data/fnf_sexDE_pval01.csv")

venn.diagram(list(ctl_desex_pval01$gene_id, fnf_desex_pval01$gene_id), 
             filename = "plots/sex_de_venn.tiff",
             category.names = c("Control", "FN-f"),
             output = TRUE,
             lwd = 0, 
             fill = c("#EDB361", "#6996C7"),
             label.col = "grey35",
             cat.col = c("#EDB361", "#6996C7"),
             cat.pos = c(10, 340), 
             cex = 2,
             cat.cex = 2.1,
             rotation.degree = 180)



both_desex_pval01 <- intersect(ctl_desex_pval01$gene_id,
                               fnf_desex_pval01$gene_id)



ctl_desex_pval01_both <- ctl_desex_pval01 %>%
  filter(gene_id %in% both_desex_pval01) %>%
  dplyr::select(gene_id, symbol, padj, log2FoldChange)
  

fnf_desex_pval01_both <- fnf_desex_pval01 %>%
  filter(gene_id %in% both_desex_pval01) %>%
  dplyr::select(gene_id, symbol, padj, log2FoldChange)

# Combine and make tidy for ggplot
ctl_fnf_desex_pval01_both <- full_join(ctl_desex_pval01_both, fnf_desex_pval01_both,
                                       by = c("gene_id", "symbol"),
                                       suffix = c("_ctl", "_fnf")) %>%
  pivot_longer(cols = c(padj_ctl, log2FoldChange_ctl, padj_fnf, log2FoldChange_fnf),
               names_to = c("stat", "condition"),
               names_sep = c("_")) %>%
  pivot_wider(names_from = "stat") %>%
  arrange(desc(abs(log2FoldChange)))


ctl_fnf_desex_pval01_both$symbol <- factor(ctl_fnf_desex_pval01_both$symbol,
                                            levels = unique(ctl_fnf_desex_pval01_both$symbol))

ggplot(ctl_fnf_desex_pval01_both, aes(x = symbol, y = log2FoldChange, fill = condition)) +
  geom_bar(stat = "identity", width = 0.7, position = "dodge") +
  scale_fill_manual(values = c("#EDB361", "#6996C7")) +
  scale_y_continuous(limits = c(-2, 12),
                     breaks = seq(-2, 12, 2),
                     name = "log~2~ Fold Change") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        panel.background = element_blank(),
        axis.line.y = element_line(color = "grey25", 
                                   linewidth = 0.25),
        axis.ticks.length.y = unit(-0.1, "cm"),
        axis.ticks = element_line(color = "grey25",
                                    linewidth = 0.25),
        legend.position = "None",
        text = element_text(family = "Helvetica"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        axis.title.y = element_markdown(),
        plot.title = ggtext::element_markdown(size = 12))  +
  labs(title = 'Effect size of sex-specific differential genes found in both
       <span style = "color:#EDB361;">**control**</span> and
       <span style = "color:#6996C7;">**FN-f**</span>')

ggsave(filename = "plots/sex_de_pval01_overlapgenes_effectsizes.pdf", 
       width = 10, height = 6, units = "in")


# untreated_sex_specific_sig <- read_csv("data/sig_untreatedsexGenes_pval01.csv") %>%
#   mutate(padj = ifelse(padj == 0, .Machine$double.xmin, padj)) %>%
#   mutate(id = "Control")
# 
# treated_sex_specific_sig <- read_csv("data/sig_treatedsexGenes_pval01.csv") %>%
#   mutate(padj = ifelse(padj == 0, .Machine$double.xmin, padj)) %>%
#   mutate(id = "FN-f")
# 
# combined_sex_specific_sig <- bind_rows(untreated_sex_specific_sig, treated_sex_specific_sig) %>%
#   group_by(gene_id) %>%
#   mutate(combined_id = paste0(id, collapse = "_")) %>%
#   ungroup()
# combined_sex_specific_sig$id <- factor(combined_sex_specific_sig$id, levels = c("Control",
#                                                                                 "FN-f"))
# 
# 
# combined_sex_specific_sig$chrom <- factor(combined_sex_specific_sig$seqnames,
#                                            levels = c(as.character(1:22), "X", "Y"))
# 
# 
# both_sex_specific_sig <- combined_sex_specific_sig %>% 
#   filter(combined_id == "Control_FN-f") %>%
#   distinct(gene_id, .keep_all = TRUE)
# 
# # Chromosome histograms of number of sig genes ----------------------------
# 
# separated <- ggplot(combined_sex_specific_sig, mapping = aes(x = chrom)) +
#   facet_wrap(vars(id), scales = "free_y", ncol = 1) +
#   scale_y_continuous(expand = c(0,0), name = "Count", limits = c(0, 40)) +
#   scale_x_discrete(name = "Chromosome", drop = FALSE) +
#   geom_histogram(stat = "count") +
#   theme_custom_general() +
#   theme(text = element_text(family = "Helvetica"),
#         axis.text.x = element_blank(),
#         axis.title = element_blank(),
#         strip.background = element_blank(),
#         strip.text = element_text(hjust = 0, face = "bold", size = 12))
# 
# overlap <- ggplot(both_sex_specific_sig, mapping = aes(x = chrom)) +
#   facet_wrap(vars(combined_id), ncol = 1, 
#              labeller = labeller(combined_id = c("Control_FN-f" = "Both"))) +
#   scale_y_continuous(expand = c(0,0), name = "Count", limits = c(0, 40)) +
#   scale_x_discrete(name = "Chromosome", drop = FALSE) +
#   geom_histogram(stat = "count") +
#   theme_custom_general() +
#   theme(text = element_text(family = "Helvetica"),
#         strip.background = element_blank(),
#         strip.text = element_text(hjust = 0, face = "bold", size = 12),
#         axis.title.y = element_blank())
# 
# sex_specific_chromhist <- wrap_elements(separated/overlap + plot_layout(heights = c(2, 1))) +
#   labs(tag = "Count")  +
#   theme(
#     plot.tag = element_text(size = rel(1), angle = 90, family = "Helvetica"),
#     plot.tag.position = "left"
#   )
# 
# ggsave(sex_specific_chromhist, filename = "plots/sex_specific_chromhist.pdf",
#        width = 8, height = 7, units = "in")
# 
# # p-value boxplots --------------------------------------------------------
# 
# separated <- ggplot(combined_sex_specific_sig, mapping = aes(x = chrom, y = -log10(padj))) +
#   facet_wrap(vars(id), scales = "free_y", ncol = 1) +
#   geom_boxplot(lwd = 0.3) +
#   scale_y_continuous(expand = c(0,0), name = "-log10(p-value)") +
#   scale_x_discrete(name = "Chromosome", drop = FALSE) +
#   theme_custom_general() +
#   theme(text = element_text(family = "Helvetica"),
#         axis.text.x = element_blank(),
#         axis.title = element_blank(),
#         strip.background = element_blank(),
#         strip.text = element_text(hjust = 0, face = "bold", size = 12))
# 
# overlap <- ggplot(both_sex_specific_sig, mapping = aes(x = chrom, y = -log10(padj))) +
#   facet_wrap(vars(combined_id), scales = "free_y", ncol = 1, 
#              labeller = labeller(combined_id = c("Control_FN-f" = "Both"))) +
#   geom_boxplot(lwd = 0.3) +
#   scale_y_continuous(expand = c(0,0), name = "-log10(p-value)") +
#   scale_x_discrete(name = "Chromosome", drop = FALSE) +
#   theme_custom_general() +
#   theme(text = element_text(family = "Helvetica"),
#         strip.background = element_blank(),
#         strip.text = element_text(hjust = 0, face = "bold", size = 12),
#         axis.title.y = element_blank())
# 
# sex_specific_pvalboxplots <- wrap_elements(separated/overlap + plot_layout(heights = c(2, 1))) +
#   labs(tag = "-log10(p-value)")  +
#   theme(
#     plot.tag = element_text(size = rel(1), angle = 90, family = "Helvetica"),
#     plot.tag.position = "left"
#   )
# 
# ggsave(sex_specific_pvalboxplots, filename = "plots/sex_specific_pvalboxplots.pdf",
#        width = 8, height = 7, units = "in")
# 
# 
# # Split by male/female and labeling with most significant gene ----------------
# 
# combined_sex_specific_sig <- combined_sex_specific_sig %>%
#   mutate(log2FC_dir = ifelse(log2FoldChange < 0, "-", "+")) %>%
#   group_by(id, chrom) %>%
#   mutate(min_gene = symbol[which.min(padj)]) %>%
#   ungroup() %>%
#   group_by(id, chrom, log2FC_dir) %>%
#   mutate(group_count = n()) %>%
#   ungroup() %>%
#   group_by(id, chrom) %>%
#   mutate(max_count = max(group_count))
# 
# 
# 
# both_sex_specific_sig <- both_sex_specific_sig %>%
#   mutate(log2FC_dir = ifelse(log2FoldChange < 0, "-", "+")) %>%
#   group_by(chrom) %>%
#   mutate(min_gene = symbol[which.min(padj)]) %>% 
#   ungroup() %>%
#   group_by(chrom, log2FC_dir) %>%
#   mutate(group_count = n()) %>%
#   ungroup() %>%
#   group_by(chrom) %>%
#   mutate(max_count = max(group_count))
#   
#   
# 
# separated <- ggplot(combined_sex_specific_sig, aes(x = chrom)) +
#   facet_wrap(vars(id), scales = "free_y", ncol = 1) +
#   geom_bar(stat = "count", position = "dodge", aes(fill = log2FC_dir)) +
#   geom_text(aes(label = min_gene, y = max_count), 
#             position = position_dodge(0.5),
#             vjust = -0.5,
#             check_overlap = TRUE,
#             size = 2.5,
#             family = "Helvetica") +
#   scale_fill_manual(values = c("#C06E8B", "#6F8CC7")) +
#   scale_x_discrete(name = "Chromosome", drop = FALSE) +
#   scale_y_continuous(expand = c(0,0), name = "Count", limits = c(0, 25)) +
#   theme_custom_general() +
#   theme(text = element_text(family = "Helvetica"),
#         legend.position = "none",
#         axis.text.x = element_blank(),
#         axis.title = element_blank(),
#         strip.background = element_blank(),
#         strip.text = element_text(hjust = 0, size = 12))
# 
# 
# 
# 
# overlap <- ggplot(both_sex_specific_sig, aes(x = chrom)) +
#   geom_bar(stat = "count", position = "dodge", aes(fill = log2FC_dir)) +
#   facet_wrap(vars(combined_id), scales = "free_y", ncol = 1, 
#              labeller = labeller(combined_id = c("Control_FN-f" = "Both"))) +
#   geom_text(aes(label = min_gene, y = max_count), 
#             position = position_dodge(0.5),
#             check_overlap = TRUE,
#             vjust = -0.5,
#             size = 2.5,
#             family = "Helvetica") +
#   scale_fill_manual(values = c("#C06E8B", "#6F8CC7")) +
#   scale_x_discrete(name = "Chromosome", drop = FALSE) +
#   scale_y_continuous(expand = c(0,0), name = "Count", limits = c(0, 25)) +
#   theme_custom_general() +
#   theme(text = element_text(family = "Helvetica"),
#         legend.position = "none",
#         strip.background = element_blank(),
#         strip.text = element_text(hjust = 0, size = 12),
#         axis.title.y = element_blank())
# 
# 
# sex_specific_siggene_chromhist <- wrap_elements(separated/overlap + plot_layout(heights = c(2, 1))) +
#   labs(tag = "Count")  +
#   plot_annotation(title = "Significant genes with sex-specific effects in control, FN-f, or both conditions",
#                   theme = theme(plot.title = element_text(family = "Helvetica"))) +
#   theme(
#     plot.tag = element_text(size = rel(1), angle = 90, family = "Helvetica"),
#     plot.tag.position = "left"
#   )
# 
# ggsave(sex_specific_siggene_chromhist, filename = "plots/sex_specific_siggene_chromhist.pdf",
#        width = 15, height = 10, units = "in")
# 
# 
# # Histogram of log2fc -----------------------------------------------------
# 
# ggplot(combined_sex_specific_sig, aes(x = abs(log2FoldChange), fill = log2FC_dir)) +
#   facet_wrap(vars(id), scales = "free_y", ncol = 1) +
#   geom_histogram() +
#   geom_vline(xintercept = 1, lty = 2, color = "grey25") +
#   scale_fill_manual(values = c("#C06E8B", "#6F8CC7")) +
#   scale_x_continuous(name = "abs(log2FC)")+
#   scale_y_continuous(expand = c(0,0)) +
#   theme_custom_general() +
#   theme(text = element_text(family = "Helvetica"),
#         legend.position = "none",
#         strip.background = element_blank(),
#         strip.text = element_text(hjust = 0, size = 12))
# 
# 
# # Volcano plot ------------------------------------------------------------
# 
# combined_sex_specific_sig <- combined_sex_specific_sig %>%
#   mutate(`-log10pval` = -log10(padj))
# 
# female_log2FClabels <- combined_sex_specific_sig %>%
#   filter(log2FC_dir == "-") %>%
#   group_by(id) %>%
#   filter(log2FoldChange == min(log2FoldChange))
#   
# female_pvallabels <- combined_sex_specific_sig %>%
#   filter(log2FC_dir == "-") %>%
#   group_by(id) %>%
#   filter(`-log10pval` == max(`-log10pval`))
#   
# male_pvallabels <- combined_sex_specific_sig %>%
#   filter(log2FC_dir == "+") %>%
#   group_by(id) %>%
#   filter(`-log10pval` != max(`-log10pval`)) %>%
#   slice_max(order_by = `-log10pval`, n = 1)
# 
# male_log2FClabels <- combined_sex_specific_sig %>%
#   filter(log2FC_dir == "+") %>%
#   group_by(id) %>%
#   filter(log2FoldChange == max(log2FoldChange))
# 
# 
#   
# filtered_labels <- combined_sex_specific_sig %>%
#   anti_join(female_log2FClabels) %>%
#   anti_join(female_pvallabels) %>%
#   anti_join(male_log2FClabels) %>%
#   anti_join(male_pvallabels)
#   
#   
#   
# 
# ggplot(combined_sex_specific_sig, aes(x = log2FoldChange, 
#                              y = `-log10pval`, 
#                              color = log2FC_dir)) +
#   facet_wrap(vars(id), scales = "free_y", ncol = 1) +
#   geom_point(size = 0.75) +
#   geom_text(data = female_log2FClabels, mapping = aes(label = symbol, 
#                                                 family = "Helvetica"),
#             size = 3,
#             nudge_x = -0.5,
#             hjust = 1) + 
#   geom_text(data = female_pvallabels, mapping = aes(label = symbol,
#                                                     family = "Helvetica"),
#             size = 3,
#             nudge_y = 5,
#             vjust = 0) +
#   geom_text(data = male_log2FClabels, mapping = aes(label = symbol,
#                                                     family = "Helvetica"),
#             size = 3,
#             nudge_x = 0.5,
#             hjust = 0) +
#   geom_text(data = male_pvallabels, mapping = aes(label = symbol,
#                                                   family = "Helvetica"),
#             size = 3,
#             nudge_y = 5,
#             vjust = 0) +
#   geom_text_repel(data = filtered_labels, aes(label = symbol, family = "Helvetica",
#                                               segment.color = "transparent"),
#                   box.padding = 0.5,
#                   size = 3,
#                   max.overlaps = 5) +
#   scale_x_continuous(breaks = seq(-10, 30, 5), limits = c(-10, 30)) +
#   scale_y_continuous(breaks = seq(0, 310, 50)) + 
#   scale_color_manual(values = c("#C06E8B", "#6F8CC7")) +
#   theme(legend.position = "none",
#         strip.background = element_blank(),
#         strip.text = element_text(hjust = 0, size = 11),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.border = element_rect(fill = "transparent", color = "grey25", linewidth = 0.5),
#         panel.background = element_blank(), 
#         axis.line = element_blank(),
#         axis.ticks.length = unit(-0.1, "cm"),
#         axis.ticks = element_line(color = "grey25",
#                                   linewidth = 0.25),
#         text = element_text(family = "Helvetica"),
#         legend.key = element_blank(),
#         legend.background = element_blank()) +
#   labs(x = "Log Fold Change",
#        y = "-log10(p-value)") 
# 
# 
# ggsave(file = "plots/sex_specific_volcano.pdf", units = "in", width = 12, height = 10)
# 
# 
# # "Venn diagram" list of genes ---------------------------------------------------
# 
# # CTL only
# ctl_only_top20_mf <- combined_sex_specific_sig %>%
#   filter(combined_id == "Control") %>%
#   group_by(log2FC_dir) %>%
#   slice_min(order_by = padj, n = 20) %>%
#   arrange(rev(padj)) %>%
#   ungroup()
# 
# # Preserve order of gene p-values
# ctl_only_top20_mf$symbol <- factor(ctl_only_top20_mf$symbol, 
#                                    levels = ctl_only_top20_mf$symbol)
# 
# 
# # FNF only  
# fnf_only_top20_mf <- combined_sex_specific_sig %>%
#   filter(combined_id == "FN-f") %>%
#   group_by(log2FC_dir) %>%
#   slice_min(order_by = padj, n = 30) %>%
#   arrange(rev(padj)) %>%
#   ungroup()
# 
# # Preserve order of gene p-values
# fnf_only_top20_mf$symbol <- factor(fnf_only_top20_mf$symbol, 
#                                    levels = fnf_only_top20_mf$symbol)
# 
# # Both only
# both_only_top20_mf <- combined_sex_specific_sig %>% 
#   filter(combined_id == "Control_FN-f") %>%
#   mutate(log2FC_dir = ifelse(log2FoldChange > 0, "+", "-")) %>%
#   group_by(id, log2FC_dir) %>%
#   slice_min(order_by = padj, n = 20) %>%
#   mutate(row = rev(row_number())) %>%
#   arrange(rev(padj)) %>%
#   ungroup() %>%
#   mutate(unique_id = row_number())
#   
# # Preserve order of gene p-values
# both_only_top20_mf$unique_id <- factor(both_only_top20_mf$unique_id, 
#                                    levels = both_only_top20_mf$unique_id)
#   
# ctl <- ggplot(ctl_only_top20_mf, aes(x = combined_id, y = symbol, fill = log2FoldChange)) +
#   geom_tile() +
#   scale_fill_steps2(low = "#C06E8B", mid = "white", high = "#6F8CC7", 
#                        limits = c(-3, 15), breaks = c(-2, -1, 0, 1, 5, 10), name = "log2FC") +
#   facet_wrap(vars(log2FC_dir), scales = "free_y", ncol = 1,
#              labeller = labeller(log2FC_dir = c("-" = "Female",
#                                                 "+" = "Male"))) +
#   guides(fill = guide_colorsteps(direction = "horizontal", title.position = "top",
#                                ticks = FALSE,
#                                title.hjust = 0.5,
#                                title.theme = element_text(size = 10, family = "Helvetica"),
#                                even.steps = TRUE, show.limits = TRUE)) +
#   geom_text(aes(label = symbol, family = "Helvetica")) +
#   theme(axis.text = element_blank(),
#         axis.title = element_blank(),
#         axis.ticks = element_blank(),
#         panel.background = element_blank(),
#         text = element_text(family = "Helvetica"),
#         strip.background = element_blank(),
#         strip.text = element_blank(),
#         legend.position = "bottom",
#         plot.title = element_text(hjust = 0.5, face = "bold")) +
#   ggtitle("Control")
# 
# fnf <- ggplot(fnf_only_top20_mf, aes(x = combined_id, y = symbol, fill = log2FoldChange)) +
#   geom_tile() +
#   scale_fill_steps2(low = "#C06E8B", mid = "white", high = "#6F8CC7", limits = c(-5, 30),
#                        breaks = c(-1, 0, 1, 10, 20), name = "log2FC") +
#   facet_wrap(vars(log2FC_dir), scales = "free_y", ncol = 1,
#              labeller = labeller(log2FC_dir = c("-" = "Female",
#                                                 "+" = "Male"))) +
#   guides(fill = guide_colorsteps(direction = "horizontal", title.position = "top",
#                                ticks = FALSE,
#                                title.hjust = 0.5,
#                                title.theme = element_text(size = 10, family = "Helvetica"),
#                                even.steps = TRUE, show.limits = TRUE)) +
#   geom_text(aes(label = symbol, family = "Helvetica")) +
#   theme(axis.text = element_blank(),
#         axis.title = element_blank(),
#         axis.ticks = element_blank(),
#         panel.background = element_blank(),
#         text = element_text(family = "Helvetica"),
#         strip.background = element_blank(),
#         strip.text = element_blank(),
#         legend.position = "bottom",
#         plot.title = element_text(hjust = 0.5, face = "bold")) +
#   ggtitle("FN-f")
# 
# both <- ggplot(both_only_top20_mf, aes(x = combined_id, y = row, fill = log2FoldChange)) +
#   geom_tile() +
#   scale_fill_steps2(low = "#C06E8B", mid = "white", high = "#6F8CC7", limits = c(-5, 10),
#                     breaks = c(-2, -1, 0, 1, 5, 10), name = "log2FC") +
#   facet_grid(log2FC_dir ~ id, scales = "free_y",
#              labeller = labeller(log2FC_dir = c("-" = "Female",
#                                                 "+" = "Male"))) +
#   guides(fill = guide_colorsteps(direction = "horizontal", title.position = "top",
#                                  ticks = FALSE,
#                                  title.hjust = 0.5,
#                                  title.theme = element_text(size = 10, family = "Helvetica"),
#                                  even.steps = TRUE, show.limits = TRUE)) +
#   geom_text(aes(label = symbol, family = "Helvetica")) +
#   theme(axis.text = element_blank(),
#         axis.title = element_blank(),
#         axis.ticks = element_blank(),
#         panel.background = element_blank(),
#         text = element_text(family = "Helvetica"),
#         strip.background = element_blank(),
#         strip.text = element_blank(),
#         legend.position = "bottom",
#         plot.title = element_text(hjust = 0.5, face = "bold")) +
#   ggtitle("Both")
# 
# ctl + both + fnf +
#   plot_layout(widths = c(0.5, 1, 0.5))
# 
# ggsave(filename = "plots/sex_specific_topGene_venndiagram.pdf", units = "in", width = 10, height = 8)
# 
# 
# 
# 
# 
# # sex-specific response to treatment l2fc heatmap ------------------------------------------------------------
# 
# load("data/dds_sex_treatmentresponse_moregenefilter.rda")
# treatmentresponse_sexDE_genes <- read_csv("data/sexDE_treatmenteffect_moregenefilter_pval05.csv") %>%
#   arrange(desc(padj))
# 
# M_l2fc <- as.data.frame(results(dds_sex_treatmentresponse_moregenefilter, 
#                                 name = "SexM.ConditionFNF")[treatmentresponse_sexDE_genes$gene_id,]) %>%
#   rownames_to_column(var = "gene_id") %>%
#   dplyr::select(gene_id, log2FoldChange) %>%
#   mutate(id = "Male")
# F_l2fc <- as.data.frame(results(dds_sex_treatmentresponse_moregenefilter, 
#                   name = "SexF.ConditionFNF")[treatmentresponse_sexDE_genes$gene_id,]) %>%
#   rownames_to_column(var = "gene_id") %>%
#   dplyr::select(gene_id, log2FoldChange) %>%
#   mutate(id = "Female")
# 
# 
# both_l2fc <- bind_rows(M_l2fc, F_l2fc)
# both_l2fc$id <- factor(both_l2fc$id, levels = c("Male", "Female"))
# both_l2fc$gene_id <- factor(both_l2fc$gene_id, levels = treatmentresponse_sexDE_genes$gene_id)
# 
# ggplot(both_l2fc, aes(x = id, y = gene_id, fill = log2FoldChange)) +
#   geom_tile() +
#   scale_fill_gradient2(low = "#73B5F9", high = "#F5D24D", mid = "black", limits = c(-6, 6)) +
#   theme(axis.title = element_blank(),
#         axis.text.y = element_blank(),
#         axis.text.x = element_text(size = 10, color = "black"),
#         axis.ticks = element_blank(),
#         panel.background = element_blank(),
#         text = element_text(family = "Helvetica"),
#         legend.title.align = 0.5,
#         legend.title = element_markdown(family = "Helvetica", size = 10),
#         plot.title = element_markdown(family = "Helvetica", hjust = 0.5, size = 10)) +
#   guides(fill = guide_colorbar(ticks.colour = NA,
#                                title = "log~2~FC")) +
#   ggtitle("Treatment log~2~FC in males and females of genes with <br> sex-specific response to FN-f (padj < 0.05)")
# ggsave(filename = "plots/sexspecifictreatment_male_female_l2fc.pdf", width = 4.5, height = 9, units = "in")
