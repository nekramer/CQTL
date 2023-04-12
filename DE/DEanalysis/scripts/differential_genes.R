library(GenomicRanges)
library(tidyverse)
library(DESeq2)
library(ggplot2)
library(ggpubr)
library(VennDiagram)
library(topGO)
library(KEGGREST)
library(ComplexHeatmap)

library(clusterProfiler)
#options(scipen=10000)

de_genes <- read_rds(file.path("data", "fnfControl_diffgenes.rds"))

# MA plot ------------------------------------------------------------
ggmaplot(as.data.frame(de_genes),
         fdr = 0.01, fc = 2, size = 0.4,
         palette = c("#EF734A", "#48A9DF", "darkgray"),
         genenames = as.vector(de_genes$symbol),
         top = 0,
         legend = "none",
         font.label = c("bold", 11),
         font.main = "bold",
         ggtheme = ggplot2::theme_minimal(),
         xlab = "log2 mean expression",
         ylab = "log2FC",
         ticks = FALSE)


ggsave("plots/2023-01-11_MAplot.pdf", units = "in",
       width = 8, height = 6)


# Comparison with Reed et al ----------------------------------------------

sig_all <- de_genes[de_genes$padj < 0.01 & abs(de_genes$log2FoldChange) > log2(1),]
reed_genes <- read_delim("data/raw/GSE150411_CRFT_differentialAnalysisFullData.txt")

differential_reed_genes <- reed_genes %>%
  filter(CLUSTER != "static") %>%
  filter(PVAL < 0.01) %>%
  filter(abs(LFC_18h) > 1)

# Filter with same p-values/log2 thresholds??? maybe based on LFC 18 h

venn.diagram(list(sig_all$gene_id, differential_reed_genes$ENSEMBL), filename = "plots/2023-01-10_Reed_overlap.tiff",
             category.names = c("Current study", "Reed et al"),
             output = TRUE,
             lwd = 0, 
             fill = c("#7bc5ee", "#fc7971"),
             label.col = "grey35",
             cat.col = NA,
             cex = 2,
             ext.pos = 90,
             ext.length = 0.6,
             margin = 0.05)


genes_both <- semi_join(as.data.frame(sig_all) %>% 
                          dplyr::rename(ENSEMBL = gene_id), 
                        differential_reed_genes, by = "ENSEMBL") 

up_both <- genes_both[which(genes_both$log2FoldChange > 0),] %>%
  # Order by decreasing log2FC
  arrange(desc(log2FoldChange))


up_both %>% dplyr::select(symbol, log2FoldChange, lfcSE, padj) %>%
  dplyr::rename(Gene = symbol) %>%
  dplyr::rename(log2FC = log2FoldChange) %>%
  dplyr::rename("FDR p-val" = padj) %>% 
  write_csv(file = "data/Reed_overlap_upgenes.csv")
  

down_both <- genes_both[which(genes_both$log2FoldChange < 0),] %>%
  # Order by decreasing log2FC
  arrange(desc(abs(log2FoldChange)))

down_both %>% dplyr::select(symbol, log2FoldChange, lfcSE, padj) %>%
  dplyr::rename(Gene = symbol) %>%
  dplyr::rename(log2FC = log2FoldChange) %>%
  dplyr::rename("FDR p-val" = padj) %>% 
  write_csv(file = "data/Reed_overlap_downgenes.csv")


genes_reed_only <- anti_join(differential_reed_genes, 
                             as.data.frame(sig_all) %>% 
                               dplyr::rename(ENSEMBL = gene_id), by = "ENSEMBL") %>%
  arrange(desc(abs(LFC_18h)))

down_late <- genes_reed_only %>% filter(CLUSTER == "down.late")
up_late <- genes_reed_only %>% filter(CLUSTER == "up.late")


genes_CQTL_only <- anti_join(as.data.frame(sig_all) %>% 
                               dplyr::rename(ENSEMBL = gene_id), 
                             differential_reed_genes, by = "ENSEMBL") %>%
  arrange(desc(abs(log2FoldChange)))

# Significant differential genes ------------------------------------------
# Define significance as padj < 0.01 abs(log2fold) change  > 2

load("data/differential_expression_dds.rda")
sig_all <- de_genes[de_genes$padj < 0.01 & abs(de_genes$log2FoldChange) > log2(2),]

write_csv(as.data.frame(sig_all), file = "data/sig_differential_genes.csv")

sig_upregulated <- de_genes[de_genes$padj < 0.01 & de_genes$log2FoldChange > log2(2),]
sig_downregulated <- sig_all[which(sig_all$log2FoldChange < 0),]

### Upregulated genes
# CXCL2
CXCL2_info <- sig_all[sig_all$symbol == "CXCL2"]
CXCL2_counts <- plotCounts(dds, gene = CXCL2_info$gene_id, intgroup = "Condition", returnData = TRUE) %>%
  rownames_to_column(var = "Sample") %>%
  separate("Sample", into = c(NA, "donor", NA, NA, NA, NA),
           sep = "_", remove = FALSE)
  

ggplot(CXCL2_counts, aes(Condition, count)) +
  geom_jitter(position = position_jitter(width = 0.1, height = 0.1)) +
  scale_y_continuous(trans = "log",
                     breaks = c(5, 20, 50, 100, 200, 1000, 5000, 10000, 50000)) +
  scale_x_discrete(labels = c("Control", "FN-f")) +
  ylab("normalized RNA-seq counts") +
  theme_light() +
  labs(title = "CXCL2", subtitle = paste0("log2FC = ", round(CXCL2_info$log2FoldChange, digits = 2))) +
  theme(axis.title.x = element_blank(), 
        plot.title = element_text(hjust = 0.5))
ggsave(filename = "plots/CXCL2_counts.pdf", units = "in",
       width = 5, height = 6)

# Donor view?
ggplot(CXCL2_counts, aes(x = donor, y = count)) +
  geom_point(aes(color = Condition)) +
  geom_line() +
  theme_minimal()

  

# IL6
IL6_info <- sig_all[sig_all$symbol == "IL6"]
IL6_counts <- plotCounts(dds, gene = IL6_info$gene_id, intgroup = "Condition", returnData = TRUE)

ggplot(IL6_counts, aes(Condition, count)) +
  geom_jitter(position = position_jitter(width = 0.1, height = 0.1)) +
  scale_y_continuous(trans = "log",
                     breaks = c(10, 100, 1000, 10000, 100000, 1000000)) +
  scale_x_discrete(labels = c("Control", "FN-f")) +
  ylab("normalized RNA-seq counts") +
  theme_light() +
  labs(title = "IL6", subtitle = paste0("log2FC = ", round(IL6_info$log2FoldChange, digits = 2))) +
  theme(axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5))

ggsave(filename = "plots/IL6_counts.pdf", units = "in",
       width = 5, height = 6)


# MMP13
MMP13_info <- sig_all[sig_all$symbol == "MMP13"]
MMP13_counts <- plotCounts(dds, gene = MMP13_info$gene_id, intgroup = "Condition", returnData = TRUE)

ggplot(MMP13_counts, aes(Condition, count)) +
  geom_jitter(position = position_jitter(width = 0.1, height = 0.1)) +
  scale_y_continuous(trans = "log",
                     breaks = c(10, 100, 1000, 10000, 100000, 1000000)) +
  scale_x_discrete(labels = c("Control", "FN-f")) +
  ylab("normalized RNA-seq counts") +
  theme_light() +
  labs(title = "MMP13", subtitle = paste0("log2FC = ", round(MMP13_info$log2FoldChange, digits = 2))) +
  theme(axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5))


ggsave(filename = "plots/MMP13_counts.pdf", units = "in",
       width = 5, height = 6)



# PLEKHM1
PLEKHM1_info <- sig_all[sig_all$symbol == "PLEKHM1"]
PLEKHM1_counts <- plotCounts(dds, gene = PLEKHM1_info$gene_id, intgroup = "Condition", returnData = TRUE) %>%
  rownames_to_column(var = "Sample") %>%
  separate("Sample", into = c(NA, "donor", NA, NA, NA, NA),
           sep = "_", remove = FALSE) 

ggplot(PLEKHM1_counts, aes(Condition, count)) +
  geom_jitter(position = position_jitter(width = 0.1, height = 0.1)) +
  scale_y_continuous(trans = "log",
                     breaks = c(0, 5, 20, 50, 200, 1000, 5000, 10000, 50000)) +
  scale_x_discrete(labels = c("Control", "FN-f")) +
  ylab("normalized RNA-seq counts") +
  theme_light() +
  labs(title = "PLEKHM1", subtitle = paste0("log2FC = ", round(PLEKHM1_info$log2FoldChange, digits = 2))) +
  theme(axis.title.x = element_blank(), 
        plot.title = element_text(hjust = 0.5))
ggsave(filename = "plots/PLEKHM1_counts.pdf", units = "in",
       width = 5, height = 6)


# WNT10B
WNT10B_info <- sig_all[sig_all$symbol == "WNT10B"]
WNT10B_counts <- plotCounts(dds, gene = WNT10B_info$gene_id, intgroup = "Condition", returnData = TRUE)

ggplot(WNT10B_counts, aes(Condition, count)) +
  geom_jitter(position = position_jitter(width = 0.1, height = 0.1)) +
  scale_y_continuous(trans = "log",
                     breaks = c(0, 1, 5, 10, 20, 50)) +
  scale_x_discrete(labels = c("Control", "FN-f")) +
  ylab("normalized RNA-seq counts") +
  theme_light() +
  labs(title = "WNT10B", subtitle = paste0("log2FC = ", round(WNT10B_info$log2FoldChange, digits = 2))) +
  theme(axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5)) 
ggsave(filename = "plots/WNT10B_counts.pdf", units = "in",
       width = 5, height = 6)


# NRP2
NRP2_info <- sig_all[sig_all$symbol == "NRP2"]
NRP2_counts <- plotCounts(dds, gene = NRP2_info$gene_id, intgroup = "Condition", returnData = TRUE)

ggplot(NRP2_counts, aes(Condition, count)) +
  geom_jitter(position = position_jitter(width = 0.1, height = 0.1)) +
  scale_y_continuous(trans = "log",
                     breaks = c(10000, 50000, 100000, 1000000)) +
  scale_x_discrete(labels = c("Control", "FN-f")) +
  ylab("normalized RNA-seq counts") +
  theme_light() +
  labs(title = "NRP2", subtitle = paste0("log2FC = ", round(NRP2_info$log2FoldChange, digits = 2))) +
  theme(axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5)) 
ggsave(filename = "plots/NRP2_counts.pdf", units = "in",
       width = 5, height = 6)

# TNKS1BP1
TNKS1BP1_info <- sig_all[sig_all$symbol == "TNKS1BP1"]
TNKS1BP1_counts <- plotCounts(dds, gene = TNKS1BP1_info$gene_id, intgroup = "Condition", returnData = TRUE)

ggplot(TNKS1BP1_counts, aes(Condition, count)) +
  geom_jitter(position = position_jitter(width = 0.1, height = 0.1)) +
  scale_y_continuous(trans = "log") +
  scale_x_discrete(labels = c("Control", "FN-f")) +
  ylab("normalized RNA-seq counts") +
  theme_light() +
  labs(title = "TNKS1BP1", subtitle = paste0("log2FC = ", round(TNKS1BP1_info$log2FoldChange, digits = 2))) +
  theme(axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5)) 
ggsave(filename = "plots/TNKS1BP1_counts.pdf", units = "in",
       width = 5, height = 6)

## Downregulated genes
# GDF5
GDF5_info <- sig_all[sig_all$symbol == "GDF5"]
GDF5_counts <- plotCounts(dds, gene = GDF5_info$gene_id, intgroup = "Condition", returnData = TRUE)

ggplot(GDF5_counts, aes(Condition, count)) +
  geom_jitter(position = position_jitter(width = 0.1, height = 0.1)) +
  scale_y_continuous(trans = "log",
                     breaks = c(50, 100, 500, 1000, 5000)) +
  scale_x_discrete(labels = c("Control", "FN-f")) +
  ylab("normalized RNA-seq counts") +
  theme_light() +
  labs(title = "GDF5", subtitle = paste0("log2FC = ", round(GDF5_info$log2FoldChange, digits = 2))) +
  theme(axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5)) 

ggsave(filename = "plots/GDF5_counts.pdf", units = "in",
       width = 5, height = 6)


# ALDH1A1
ALDH1A1_info <- sig_all[sig_all$symbol == "ALDH1A1"]
ALDH1A1_counts <- plotCounts(dds, gene = ALDH1A1_info$gene_id, intgroup = "Condition", returnData = TRUE)

ggplot(ALDH1A1_counts, aes(Condition, count)) +
  geom_jitter(position = position_jitter(width = 0.1, height = 0.1)) +
  scale_y_continuous(trans = "log",
                     breaks = c(0, 2, 5, 10, 50, 100, 250)) +
  scale_x_discrete(labels = c("Control", "FN-f")) +
  ylab("normalized RNA-seq counts") +
  theme_light() +
  labs(title = "ALDH1A1", subtitle = paste0("log2FC = ", round(ALDH1A1_info$log2FoldChange, digits = 2), ", FDR p-val = ", round(ALDH1A1_info$padj))) +
  theme(axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5)) 

ggsave(filename = "plots/ALDH1A1_counts.pdf", units = "in",
       width = 5, height = 6)


# COL2A1
COL2A1_info <- sig_all[sig_all$symbol == "COL2A1"]
COL2A1_counts <- plotCounts(dds, gene = COL2A1_info$gene_id, intgroup = "Condition", returnData = TRUE)

ggplot(COL2A1_counts, aes(Condition, count)) +
  geom_jitter(position = position_jitter(width = 0.1, height = 0.1)) +
  scale_y_continuous(trans = "log",
                     breaks = c(0, 1000, 10000, 25000, 50000, 100000, 500000)) +
  scale_x_discrete(labels = c("Control", "FN-f")) +
  ylab("normalized RNA-seq counts") +
  theme_light() +
  labs(title = "COL2A1", subtitle = paste0("log2FC = ", round(COL2A1_info$log2FoldChange, digits = 2))) +
  theme(axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5)) 

ggsave(filename = "plots/COL2A1_counts.pdf", units = "in",
       width = 5, height = 6)


# GALNT8
GALNT8_info <- sig_all[sig_all$symbol == "GALNT8"]
GALNT8_counts <- plotCounts(dds, gene = GALNT8_info$gene_id, intgroup = "Condition", returnData = TRUE)

ggplot(GALNT8_counts, aes(Condition, count)) +
  geom_jitter(position = position_jitter(width = 0.1, height = 0.1)) +
  scale_y_continuous(trans = "log") +
  scale_x_discrete(labels = c("Control", "FN-f")) +
  ylab("normalized RNA-seq counts") +
  theme_light() +
  labs(title = "GALNT8", subtitle = paste0("log2FC = ", round(GALNT8_info$log2FoldChange, digits = 2))) +
  theme(axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5)) 

# Overlap with Boer high confidence effector genes
boer_genes <- read_csv("data/raw/Boer_effector_genes.csv")
boer_overlap <- sig_all[which(sig_all$symbol %in% boer_genes$Gene),]

# GO/KEGG -----------------------------------------------------------------
geneList_up <- sig_upregulated$padj
names(geneList_up) <- sig_upregulated$symbol

GOdata <- new("topGOdata",
              ontology = "BP",
              allGenes = geneList_up,
              geneSelectionFun = function(x)x,
              annot = annFUN.org, mapping = "org.Hs.eg.db", ID = "symbol")

resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")

go_table_up <- GenTable(GOdata, raw.p.value = resultFisher, topNodes = length(resultFisher@score),
                            numChar = 120) %>%
  arrange(raw.p.value)

sig_go_up <- go_table_up %>%
  filter(raw.p.value < 0.05)

geneList_down <- sig_downregulated$padj
names(geneList_down) <- sig_downregulated$symbol

GOdata <- new("topGOdata",
              ontology = "BP",
              allGenes = geneList_down,
              geneSelectionFun = function(x)x,
              annot = annFUN.org, mapping = "org.Hs.eg.db", ID = "symbol")

resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")

go_table_down <- GenTable(GOdata, raw.p.value = resultFisher, topNodes = length(resultFisher@score),
                        numChar = 120) %>%
  arrange(raw.p.value)

sig_go_down <- go_table_down %>%
  filter(raw.p.value < 0.05)



keggPathways <- keggList("pathway", "hsa")
geneList <- sig_all$padj
names(geneList) <- sig_all$symbol


# Pull all genes for each pathway
pathway.codes <- sub("path:", "", names(keggPathways)) 
genes.by.pathway <- lapply(pathway.codes,
                           function(pwid){
                             print(pwid)
                             pw <- keggGet(pwid)
                             if (is.null(pw[[1]]$GENE)) return(NA)
                             pw2 <- pw[[1]]$GENE[c(FALSE,TRUE)]
                             pw2 <- unlist(lapply(strsplit(pw2, split = ";", fixed = T), function(x)x[1]))
                             return(pw2)
                           }
)
names(genes.by.pathway) <- pathway.codes



pVals.by.pathway <- t(sapply(names(genes.by.pathway),
                             function(pathway) {
                               pathway.genes <- genes.by.pathway[[pathway]]
                               list.genes.in.pathway <- intersect(names(geneList), pathway.genes)
                               list.genes.not.in.pathway <- setdiff(names(geneList), list.genes.in.pathway)
                               scores.in.pathway <- geneList[list.genes.in.pathway]
                               scores.not.in.pathway <- geneList[list.genes.not.in.pathway]
                               if (length(scores.in.pathway) > 0){
                                 p.value <- wilcox.test(scores.in.pathway, scores.not.in.pathway, alternative = "less")$p.value
                               } else{
                                 p.value <- NA
                               }
                               return(c(p.value = p.value, Annotated = length(list.genes.in.pathway)))
                             }
))

outdat <- data.frame(pathway.code = rownames(pVals.by.pathway))
outdat$pathway.name <- keggPathways[paste0("path:", outdat$pathway.code)]
outdat$p.value <- pVals.by.pathway[,"p.value"]
outdat$Annotated <- pVals.by.pathway[,"Annotated"]
outdat <- outdat[order(outdat$p.value),]
outdat$pathway.name <- gsub("- Homo sapiens \\(human\\)","", outdat$pathway.name)

kegg_plotting <- outdat %>%
  slice_min(order_by = p.value, n = 20) %>%
  mutate(logpval = -1* log10(p.value)) %>%
  arrange(p.value)

kegg_plotting$pathway.name <- factor(kegg_plotting$pathway.name, levels = kegg_plotting$pathway.name)

ggplot(kegg_plotting, aes(x = 1, y = pathway.name, fill = logpval)) +
  geom_tile() +
  scale_y_discrete(limits = rev) +
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank()) +
  labs(fill = "-log10 p-value")
ggsave("plots/keggpathways.pdf", units = "in",
       width = 5.5,
       height = 5)


# Heatmap -----------------------------------------------------------------

de_genes <- read_rds(file.path("data", "fnfControl_diffgenes.rds"))

load("data/differential_expression_dds.rda")
#sig_all <- de_genes[de_genes$padj < 0.01 & abs(de_genes$log2FoldChange) > log2(2),]

normCounts <- counts(dds, normalized = TRUE)[de_genes$gene_id,]

Heatmap(normCounts, cluster_rows = TRUE, cluster_columns = TRUE)

