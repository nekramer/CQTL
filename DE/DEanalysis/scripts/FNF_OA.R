library(DESeq2)
library(tidyverse)
library(ggtext)
library(colorspace)
source("scripts/plotting_utils.R")

fit_de_model <- function(dds, filter = 0.05){
  ## Filter out lowly expressed genes
  # 10 reads in at least filter% of samples
  keep <- rowSums(counts(dds) >= 10) >= ceiling(nrow(colData(dds))*filter)
  dds <- dds[keep,]

  ## Fit model
  dds_fit <- DESeq(dds)
  return(dds_fit)
}

### OA and FNF together -----------------------------------------------------


load("data/2023-07-19_gse_oa_fnf.rda")

dds_oa_fnf <- DESeqDataSet(gse_oa_fnf, design = ~Sex + Age + Race + Condition)
dds_oa_fnf <- fit_de_model(dds_oa_fnf)
save(dds_oa_fnf, file = "data/dds_oa_fnf.rda")

# PCA
# Normalized gene counts
normCounts <- t(assay(vst(dds_oa_fnf)))

normCounts <- as.data.frame(normCounts) %>%
  filter(!grepl("FNF", rownames(normCounts)))
# Calculate principal components and percent variance explained
pcs_normCounts <- prcomp(normCounts)
percentVar <- pcs_normCounts$sdev^2/sum(pcs_normCounts$sdev^2)


# Create dataframe for plotting
pc_df <- data.frame("PC1" = pcs_normCounts$x[,1],
                    "PC2" = pcs_normCounts$x[,2]) %>%
  rownames_to_column(var = "Sample") %>%
  mutate(group = ifelse(grepl("OA", Sample), "OA", gsub("^.*_(CTL|FNF)_.*$", "\\1", Sample)))

ggplot(pc_df, aes(x = PC1, y = PC2, color = group)) +
  geom_point()



# Gene direction comparison -----------------------------------------------

oa_genes_shrink <- lfcShrink(dds_oa_fnf, coef = "Condition_OA_vs_CTL") %>%
  as.data.frame() %>%
  rownames_to_column(var = "gene_id")

oa_genes_shrink <- inner_join(x = oa_genes_shrink,
                              y = as.data.frame(rowData(gse_oa_fnf)) %>%
                                dplyr::select(c("gene_id", "symbol", "tx_ids")),
                              by = "gene_id")
write_csv(oa_genes_shrink, file  = "data/oa_genes_shrink.csv")  


## P-val < 0.01 and abs(l2fc) > 1
oa_genes_sig_pval01_l2fc1 <- oa_genes_shrink %>%
  filter(padj < 0.01 & abs(log2FoldChange) > 1)
write_csv(oa_genes_sig_pval01_l2fc1, file = "data/sig_oaGenes_pval01_l2fc1.csv")
oa_genes_sig_pval01_l2fc1_up <- oa_genes_sig_pval01_l2fc1 %>%
  filter(log2FoldChange > 0) %>%
  arrange(desc(log2FoldChange))
oa_genes_sig_pval01_l2fc1_down <- oa_genes_sig_pval01_l2fc1 %>%
  filter(log2FoldChange < 0) %>%
  arrange(log2FoldChange)

## P-val < 0.01 and abs(l2fc) > 1, top 100 l2fc
oa_genes_sig_pval01_l2fc1_up_100 <- oa_genes_sig_pval01_l2fc1_up %>%
  slice_head(n = 100)
oa_genes_sig_pval01_l2fc1_down_100 <- oa_genes_sig_pval01_l2fc1_down %>%
  slice_head(n = 100)



## P-val < 0.01 and abs(l2fc) > 2
oa_genes_sig_pval01_l2fc2 <- oa_genes_shrink %>%
  filter(padj < 0.01 & abs(log2FoldChange) > 2)
oa_genes_sig_pval01_l2fc2_up <- oa_genes_sig_pval01_l2fc2 %>%
  filter(log2FoldChange > 0)
oa_genes_sig_pval01_l2fc2_down <- oa_genes_sig_pval01_l2fc2 %>%
  filter(log2FoldChange < 0)

#fnf_genes <- read_csv("data/de_genes_results.csv")
fnf_genes_shrink <- lfcShrink(dds_oa_fnf, coef = "Condition_FNF_vs_CTL") %>%
  as.data.frame() %>%
  rownames_to_column(var = "gene_id")

fnf_genes_shrink <- inner_join(x = fnf_genes_shrink,
                              y = as.data.frame(rowData(gse_oa_fnf)) %>%
                                dplyr::select(c("gene_id", "symbol", "tx_ids")),
                              by = "gene_id")

write_csv(fnf_genes_shrink, file = "data/fnf_genes_shrink.csv")

fnf_genes_oa_up_subset <- fnf_genes_shrink %>%
  filter(gene_id %in% oa_genes_sig_pval01_l2fc1_up$gene_id) %>%
  mutate(group = "Up in OA")

fnf_genes_oa_down_subset <- fnf_genes_shrink %>%
  filter(gene_id %in% oa_genes_sig_pval01_l2fc1_down$gene_id) %>%
  mutate(group = "Down in OA")

# fnf_genes_oa_up_subset <- fnf_genes %>%
#   filter(gene_id %in% oa_genes_sig_pval01_l2fc1_up$gene_id) %>%
#   mutate(group = "Up in OA")
# 
# fnf_genes_oa_down_subset <- fnf_genes %>%
#   filter(gene_id %in% oa_genes_sig_pval01_l2fc1_down$gene_id) %>%
#   mutate(group = "Down in OA")


up_test <- wilcox.test(x = fnf_genes_oa_up_subset$log2FoldChange,
            y = fnf_genes_shrink %>% 
              filter(!gene_id %in% fnf_genes_oa_up_subset$gene_id) %>% 
              pull(log2FoldChange),
            alternative = "greater")

down_test <- wilcox.test(x = fnf_genes_oa_down_subset$log2FoldChange,
                       y = fnf_genes_shrink %>% 
                         filter(!gene_id %in% fnf_genes_oa_down_subset$gene_id) %>% 
                         pull(log2FoldChange),
                       alternative = "less")


fnf_genes_oa_subset <- bind_rows(fnf_genes_oa_up_subset,
                                 fnf_genes_oa_down_subset)
fnf_genes_oa_subset$group <- factor(fnf_genes_oa_subset$group, 
                                    levels = c("Up in OA", "Down in OA"))

ggplot(fnf_genes_oa_subset, aes(x = group, y = log2FoldChange, 
                                color = group, fill = group)) +
  scale_fill_manual(values = c(log2fcColors[["+"]], log2fcColors[["-"]])) +
  scale_color_manual(values = c(darken(log2fcColors[["+"]], 0.3), 
                                darken(log2fcColors[["-"]], 0.3))) +
  scale_y_continuous(name = "log~2~(fold change)<br>in response to FN-f") +
  geom_hline(yintercept = 0, lty = 2) +
  geom_boxplot(outlier.shape = 1, width = 0.35) +
  stat_boxplot(geom = "errorbar", width = 0.2) +
  # annotate("text",
  #          label = paste0("pval = ", format(up_test$p.value, digits = 3)),
  #          x = 1, y = Inf, family = "Helvetica", size = 3.5, vjust = 1) +
  # annotate("text",
  #          label = paste0("pval = ",format(down_test$p.value, digits = 3)),
  #          x = 2, y = Inf, family = "Helvetica", size = 3.5, vjust = 1) +
  annotate("text", label = "****", x = 2, y = Inf, family = "Helvetica",
           size = 5, vjust = 1, fontface = "bold") + 
  annotate("richtext", label = 'p-value <b>****</b> = < 2.2e-16',
           x = 0.5, y = Inf, hjust = 0, vjust = 1, family = "Helvetica",
           size = 3, label.color = NA) +
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


ggsave(file = "plots/OA_FNF_boxplots.pdf", width = 6, height = 6, units = "in")


# RAAK genes --------------------------------------------------------------

raak_genes <- read_csv("data/RAAK_genes.csv",
                       col_select = c("ENSEMBL", "HGNC", "RAAK_PVAL",
                                      "RAAK_FC", "RAAK_LFC"))

fnf_genes <- read_csv("data/de_genes_results.csv")


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
  # annotate("text", label = "****", x = 1, y = Inf, family = "Helvetica",
  #          size = 5, vjust = 1, fontface = "bold") +
  # annotate("text", label = "****", x = 2, y = Inf, family = "Helvetica",
  #          size = 5, vjust = 1, fontface = "bold") +
  # annotate("richtext", label = 'p-value <b>****</b> = < 0.0001',
  #          x = 0.5, y = -7, hjust = 0, vjust = 1, family = "Helvetica",
  #          size = 3.5, label.color = NA) +
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

ggsave(file = "plots/RAAKOA_FNF_boxplots.pdf", width = 6, height = 6, units = "in")

#### OA alone (with Control) -------------------------------------------------
# 
# load("data/2023-07-20_gse_oa.rda")
# 
# dds_oa <- DESeqDataSet(gse_oa, design = ~Sex + Age + Race + Condition)
# dds_oa <- fit_de_model(dds_oa)
# save(dds_oa, file = "data/dds_oa.rda")
# 
# oa_shrink <- lfcShrink(dds_oa, coef = "Condition_OA_vs_CTL") %>%
#   as.data.frame() %>%
#   rownames_to_column(var = "gene_id")
# 
# oa_shrink <- inner_join(x = oa_shrink,
#                               y = as.data.frame(rowData(gse_oa)) %>%
#                                 dplyr::select(c("gene_id", "symbol", "tx_ids")),
#                               by = "gene_id")
# write_csv(oa_shrink, file  = "data/oa_shrink.csv")
# 
# 
# 
# oa_sig_pval01_l2fc1 <- oa_shrink %>%
#   filter(padj < 0.01 & abs(log2FoldChange) > 1)
# 
# oa_sig_pval01_l2fc1_up <- oa_sig_pval01_l2fc1 %>%
#   filter(log2FoldChange > 0)
# 
# oa_sig_pval01_l2fc1_down <- oa_sig_pval01_l2fc1 %>%
#   filter(log2FoldChange < 0)
# 
# 
# fnf_genes <- read_csv("data/de_genes_results.csv")
# fnf_oa_up_subset <- fnf_genes %>%
#   filter(gene_id %in% oa_sig_pval01_l2fc1_up$gene_id) %>%
#   mutate(group = "Up in OA")
# 
# fnf_oa_down_subset <- fnf_genes %>%
#   filter(gene_id %in% oa_sig_pval01_l2fc1_down$gene_id) %>%
#   mutate(group = "Down in OA")
# 
# 
# up_test2 <- wilcox.test(x = oa_sig_pval01_l2fc1_up$log2FoldChange,
#                        y = fnf_oa_up_subset$log2FoldChange,
#                        alternative = "two.sided")
# 
# down_test2 <- wilcox.test(x = oa_sig_pval01_l2fc1_down$log2FoldChange,
#                          y = fnf_oa_down_subset$log2FoldChange,
#                          alternative = "two.sided")
# 
# 
# fnf_oa_subset <- bind_rows(fnf_oa_up_subset,
#                                  fnf_oa_down_subset)
# fnf_oa_subset$group <- factor(fnf_oa_subset$group, 
#                                     levels = c("Up in OA", "Down in OA"))
# 
# ggplot(fnf_oa_subset, aes(x = group, y = log2FoldChange, 
#                                 color = group, fill = group)) +
#   scale_fill_manual(values = c(log2fcColors[["+"]], log2fcColors[["-"]])) +
#   geom_hline(yintercept = 0, lty = 2) +
#   geom_boxplot() +
#   theme(panel.background = element_blank(),
#         axis.ticks.x = element_blank(),
#         axis.title.x = element_blank(),
#         legend.position = "None")
