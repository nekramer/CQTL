library(ggrepel)
library(ggtext)
source("scripts/plotting_utils.R")

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


fnf_de_genes <- read_csv("data/sig_deGenes_pval01_l2fc1.csv")

oa_de_genes <- read_csv("data/RAAK_genes.csv",
                        col_select = c("ENSEMBL", "HGNC", "RAAK_PVAL",
                                       "RAAK_FC", "RAAK_LFC")) %>%
  dplyr::rename(symbol = HGNC) %>%
  dplyr::rename(log2FoldChange = RAAK_LFC) %>%
  dplyr::rename(gene_id = ENSEMBL)



fnf_oa_de_genes <- fnf_de_genes %>%
  filter(symbol %in% oa_de_genes$symbol)


# Get age de genes only differential in FN-f and get l2fc in FN-f
age_fnf_degenes <- union_sig_genes %>%
  filter(gene_id %in% fnf_de_genes$gene_id & !gene_id %in% oa_de_genes$gene_id) %>%
  left_join(fnf_de_genes %>% dplyr::select(gene_id, symbol, log2FoldChange), by = "gene_id") %>%
  dplyr::rename(fnf_log2FoldChange = log2FoldChange)

# Get age de genes only differential in OA (and get OA l2fc)
age_oa_degenes <- union_sig_genes %>%
  filter(!gene_id %in% fnf_de_genes$gene_id & gene_id %in% oa_de_genes$gene_id) %>%
  left_join(oa_de_genes %>% dplyr::select(gene_id, symbol, log2FoldChange), by = "gene_id") %>%
  dplyr::rename(oa_log2FoldChange = log2FoldChange)

# Age DE genes differential in both
age_fnfoa_degenes <- union_sig_genes %>%
  filter(gene_id %in% fnf_de_genes$gene_id & gene_id %in% oa_de_genes$gene_id) %>%
  left_join(fnf_de_genes %>% dplyr::select(gene_id, symbol, log2FoldChange), by = "gene_id") %>%
  dplyr::rename(fnf_log2FoldChange = log2FoldChange) %>%
  left_join(oa_de_genes %>% dplyr::select(gene_id, log2FoldChange), by = "gene_id") %>%
  dplyr::rename(oa_log2FoldChange = log2FoldChange)



no_group <- union_sig_genes %>%
  filter(!gene_id %in% fnf_de_genes$gene_id & !gene_id %in% oa_de_genes$gene_id)



# Calculate enrichment p-value with permutation test ----------------------
de_genes_results <- read_csv("data/de_genes_results.csv")

nPerm <- 1000
# FN-f group
fnf_randomOverlaps <- c() 

for (i in 1:nPerm){
  # Randomly select number of sex_degenes_union from all genes
  random_genes <- de_genes_results[sample(nrow(de_genes_results),
                                          nrow(union_sig_genes)),]
  
  # Check how many overlap with FN-f
  fnfOverlap <- length(which(random_genes$symbol %in% fnf_de_genes$symbol))
  
  
  fnf_randomOverlaps[i] <- fnfOverlap
  
}
# Proportion of fnf_randomOverlaps that are greater than sex_fnf_degenes overlap
fnf_pval <- length(which(fnf_randomOverlaps > nrow(age_fnf_degenes)))/nPerm



oa_randomOverlaps <- c()

for (i in 1:nPerm){
  
  # Randomly select number of sex_degenes_union from all genes
  random_genes <- de_genes_results[sample(nrow(de_genes_results),
                                          nrow(union_sig_genes)),]
  
  # Check how many overlap with RAAK
  oaOverlap <- length(which(random_genes$symbol %in% oa_de_genes$symbol))
  
  oa_randomOverlaps[i] <- oaOverlap
  
}

oa_pval <- length(which(oa_randomOverlaps > nrow(age_oa_degenes)))/nPerm 

fnf_oa_randomOverlaps <- c()
for (i in 1:nPerm){
  
  # Randomly select number of sex_degenes_union from all genes
  random_genes <- de_genes_results[sample(nrow(de_genes_results),
                                          nrow(union_sig_genes)),]
  
  # Check how many overlap with both FN-f and RAAK
  fnf_oaOverlap <- length(which(random_genes$symbol %in% fnf_oa_de_genes$symbol))
  
  
  fnf_oa_randomOverlaps[i] <- fnf_oaOverlap
  
}
fnf_oa_pval <- length(which(fnf_oa_randomOverlaps >= nrow(age_fnfoa_degenes)))/nPerm


# Pie chart ---------------------------------------------------------------


pie_chart_data <- data.frame(group = c("FN-f", "FN-f and OA", "OA", "NA"),
                             number = c(nrow(age_fnf_degenes),
                                        nrow(age_fnfoa_degenes),
                                        nrow(age_oa_degenes),
                                        nrow(no_group)),
                             pval = c(fnf_pval, fnf_oa_pval, oa_pval, NA))
pie_chart_data$group <- factor(pie_chart_data$group,
                               levels = c("OA", "FN-f and OA", "FN-f", "NA"))
ggplot(pie_chart_data, aes(y = number, x = "", fill = group)) +
  geom_bar(stat = "identity", width = 1) +
  # Actual numbers
  geom_text(aes(x = 1.3, label = ifelse(group == "NA", NA, number)), 
            position = position_stack(vjust = 0.5),
            family = "Helvetica",
            fontface = "bold",
            size = 5) +
  # Category names
  geom_text(aes(x = c(1.68, 1.56, 1.57, NA),
                y = c(70, 76, 78.5, NA),
                label = group,
                color = group), hjust = 0,
            family = "Helvetica",
            fontface = "bold",
            size = 6) +
  # p values
  geom_text(aes(x = c(1.65, 1.56, 1.58, NA),
                y = c(70.5, 76.5, 79, NA),
                label = ifelse(group == "NA",
                               NA,
                               paste0("pval = ", pval)),
                color = group), hjust = 0,
            vjust = 1,
            family = "Helvetica",
            size = 4) +
  scale_fill_manual(values = c(
    "FN-f and OA" = "#47B39C", 
    "FN-f"= "#FFC154", 
    "NA" = "grey85", 
    "OA" = "#2D87BB")) +
  scale_color_manual(values = c(
    "FN-f and OA" = "#47B39C", 
    "FN-f"= "#FFC154", 
    "NA" = "grey85", 
    "OA" = "#2D87BB")) +
  coord_polar("y", start = 6.9*pi/12, clip = 'off') +
  theme_void() +
  theme(legend.position = "None",
        plot.title = element_text(family = "Helvetica", hjust = 0.5, vjust = -7,
                                  face = "bold", size = 18)) +
  ggtitle("Overlap with differential OA and FN-f genes")

ggsave(filename = "plots/age_FNF_OA_piechart.pdf",
       width = 8, height = 6.5, unit = "in",
       bg = "transparent")
