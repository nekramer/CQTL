library(tidyverse)
library(ggrepel)
library(ggtext)
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


# Sex gene overlaps -------------------------------------------------------

ctl_sex_degenes <- read_csv("data/ctl_sexDE_pval01.csv")
fnf_sex_degenes <- read_csv("data/fnf_sexDE_pval01.csv")

sex_degenes_union <- full_join(ctl_sex_degenes %>%
            dplyr::select("symbol", "gene_id", "log2FoldChange") %>%
              mutate(condition = "ctl"),
          fnf_sex_degenes %>%
            dplyr::select("symbol", "gene_id", "log2FoldChange") %>%
            mutate(condition = "fnf")) %>%
  distinct(gene_id, .keep_all = TRUE) %>%
  dplyr::rename(sex_log2FoldChange = log2FoldChange)
  # group_by(symbol, gene_id) %>%
  # summarise(adj_log2FoldChange = mean(log2FoldChange))



fnf_de_genes <- read_csv("data/sig_deGenes_pval01_l2fc1.csv")

oa_de_genes <- read_csv("data/RAAK_genes.csv",
                       col_select = c("ENSEMBL", "HGNC", "RAAK_PVAL",
                                      "RAAK_FC", "RAAK_LFC")) %>%
  dplyr::rename(symbol = HGNC) %>%
  dplyr::rename(log2FoldChange = RAAK_LFC)


fnf_oa_de_genes <- fnf_de_genes %>%
  filter(symbol %in% oa_de_genes$symbol)


# Get sex de genes only differential in FN-f (and get FN-f l2fc)
sex_fnf_degenes <- sex_degenes_union %>%
  filter(symbol %in% fnf_de_genes$symbol & !symbol %in% oa_de_genes$symbol) %>%
  left_join(fnf_de_genes %>% dplyr::select(symbol, log2FoldChange)) %>%
  dplyr::rename(fnf_log2FoldChange = log2FoldChange) %>%
  mutate(sex_group = ifelse(sex_log2FoldChange < 0, "female", "male")) %>%
  mutate(fnf_group = ifelse(fnf_log2FoldChange < 0, "down", "up")) %>%
  mutate(sex_fnf_group = paste0(sex_group, "_", fnf_group))

sex_fnf_degenes_all <- sex_degenes_union %>%
  filter(symbol %in% fnf_de_genes$symbol) %>%
  left_join(fnf_de_genes %>% dplyr::select(symbol, log2FoldChange)) %>%
  dplyr::rename(fnf_log2FoldChange = log2FoldChange) %>%
  mutate(sex_group = ifelse(sex_log2FoldChange < 0, "female", "male")) %>%
  mutate(fnf_group = ifelse(fnf_log2FoldChange < 0, "down", "up")) %>%
  mutate(sex_fnf_group = paste0(sex_group, "_", fnf_group))

# Get sex de genes only differential in OA (and get OA l2fc)
sex_oa_degenes <- sex_degenes_union %>%
  filter(!symbol %in% fnf_de_genes$symbol & symbol %in% oa_de_genes$symbol) %>%
  left_join(oa_de_genes %>% dplyr::select(symbol, log2FoldChange)) %>%
  dplyr::rename(oa_log2FoldChange = log2FoldChange) %>%
  mutate(sex_group = ifelse(sex_log2FoldChange < 0, "female", "male")) %>%
  mutate(oa_group = ifelse(oa_log2FoldChange < 0, "down", "up")) %>%
  mutate(sex_oa_group = paste0(sex_group, "_", oa_group))

sex_oa_degenes_all <- sex_degenes_union %>%
  filter(symbol %in% oa_de_genes$symbol) %>%
  left_join(oa_de_genes %>% dplyr::select(symbol, log2FoldChange)) %>%
  dplyr::rename(oa_log2FoldChange = log2FoldChange) %>%
  mutate(sex_group = ifelse(sex_log2FoldChange < 0, "female", "male")) %>%
  mutate(oa_group = ifelse(oa_log2FoldChange < 0, "down", "up")) %>%
  mutate(sex_oa_group = paste0(sex_group, "_", oa_group))

sex_fnfoa_degenes <- sex_degenes_union %>%
  filter(symbol %in% fnf_de_genes$symbol & symbol %in% oa_de_genes$symbol) %>%
  left_join(fnf_de_genes %>% dplyr::select(symbol, log2FoldChange)) %>%
  dplyr::rename(fnf_log2FoldChange = log2FoldChange) %>%
  left_join(oa_de_genes %>% dplyr::select(symbol, log2FoldChange)) %>%
  dplyr::rename(oa_log2FoldChange = log2FoldChange) %>%
  mutate(sex_group = ifelse(sex_log2FoldChange < 0, "female", "male")) %>%
  mutate(fnf_group = ifelse(fnf_log2FoldChange < 0, "down", "up")) %>%
  mutate(oa_group = ifelse(oa_log2FoldChange < 0, "down", "up"))



no_group <- sex_degenes_union %>%
  filter(!symbol %in% fnf_de_genes$symbol & !symbol %in% oa_de_genes$symbol)


# Calculate enrichment p-value with permutation test ----------------------
de_genes_results <- read_csv("data/de_genes_results.csv")

# FN-f group
fnf_randomOverlaps <- c() 
for (i in 1:nPerm){
  # Randomly select number of sex_degenes_union from all genes
  random_genes <- de_genes_results[sample(nrow(de_genes_results),
                                          nrow(sex_degenes_union)),]
  
  # Check how many overlap with FN-f
  
  fnfOverlap <- length(which(random_genes$symbol %in% fnf_de_genes$symbol))
  
  
  fnf_randomOverlaps[i] <- fnfOverlap
  
}
# Proportion of fnf_randomOverlaps that are greater than sex_fnf_degenes overlap

fnf_pval <- length(which(fnf_randomOverlaps > nrow(sex_fnf_degenes)))/nPerm

oa_randomOverlaps <- c()
for (i in 1:nPerm){
  
  # Randomly select number of sex_degenes_union from all genes
  random_genes <- de_genes_results[sample(nrow(de_genes_results),
                                          nrow(sex_degenes_union)),]
  
  # Check how many overlap with RAAK
  oaOverlap <- length(which(random_genes$symbol %in% oa_de_genes$symbol))
  
  oa_randomOverlaps[i] <- oaOverlap
  
}

oa_pval <- length(which(oa_randomOverlaps > nrow(sex_oa_degenes)))/nPerm 

fnf_oa_randomOverlaps <- c()
for (i in 1:nPerm){
  
  # Randomly select number of sex_degenes_union from all genes
  random_genes <- de_genes_results[sample(nrow(de_genes_results),
                                          nrow(sex_degenes_union)),]
  
  # Check how many overlap with both FN-f and RAAK
  fnf_oaOverlap <- length(which(random_genes$symbol %in% fnf_oa_de_genes$symbol))
  
  
  fnf_oa_randomOverlaps[i] <- fnf_oaOverlap

}

fnf_oa_pval <- length(which(fnf_oa_randomOverlaps >= nrow(sex_fnfoa_degenes)))/nPerm


overlapData <- data.frame(group = c("FN-f", "FN-f and OA", "OA", "NA"),
                          number = c(nrow(sex_fnf_degenes),
                                     nrow(sex_fnfoa_degenes),
                                     nrow(sex_oa_degenes),
                                     nrow(no_group)),
                          pval = c(fnf_pval, fnf_oa_pval, oa_pval, NA),
                          expected_value = c(median(fnf_randomOverlaps),
                                             median(fnf_oa_randomOverlaps),
                                             median(oa_randomOverlaps),
                                             NA),
                          expected_value_q1 = c(quantile(fnf_randomOverlaps, probs = 0.25),
                                                quantile(fnf_oa_randomOverlaps, probs = 0.25),
                                                quantile(oa_randomOverlaps, probs = 0.25),
                                                NA),
                          expected_value_q3 = c(quantile(fnf_randomOverlaps, probs = 0.75),
                                                quantile(fnf_oa_randomOverlaps, probs = 0.75),
                                                quantile(oa_randomOverlaps, probs = 0.75),
                                                NA))
# Pie chart ---------------------------------------------------------------

pie_chart_data <- overlapData
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
                y = c(95, 102.75, 105.9, NA),
                label = group,
                color = group), hjust = 0, 
            family = "Helvetica",
            fontface = "bold",
            size = 6) +
  # p values
  geom_text(aes(x = c(1.65, 1.56, 1.58, NA),
                y = c(95.75, 103.75, 106.75, NA),
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

ggsave(filename = "plots/sex_FNF_OA_piechart.pdf",
       width = 8, height = 6.25, unit = "in",
       bg = "transparent")


# Bar plot of overlaps ----------------------------------------------------

barplotData <- overlapData %>% 
  filter(group != "NA") %>%
  dplyr::rename(observed = number) %>%
  pivot_longer(cols = c("observed", "expected_value"),
               names_to = "observed_expected",
               values_to = "num") %>%
  mutate(expected_value_q1 = ifelse(observed_expected == "expected_value", expected_value_q1, NA)) %>%
  mutate(expected_value_q3 = ifelse(observed_expected == "expected_value", expected_value_q3, NA)) %>%
  mutate(pval = ifelse(observed_expected == "expected_value", NA, pval))
barplotData$group <- factor(barplotData$group, levels = c("FN-f", "OA", "FN-f and OA"))


ggplot(barplotData, aes(x = group, y = num, fill = observed_expected, color = observed_expected, group = observed_expected)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("#FFC154", "#47B39C"), labels = c("Expected", "Observed")) +
  scale_color_manual(values = c("#FFC154", "#47B39C"), labels = c("Expected", "Observed")) +
  scale_y_continuous(name = "Number of overlapping sex-specific genes", limits = c(0, 30), expand = c(0,0)) +
  geom_errorbar(aes(ymin = expected_value_q1, ymax = expected_value_q3), 
                color = "black",
                position = position_dodge(width = 0.9), width = 0.25,
                linewidth = 0.25) +
  coord_cartesian(clip = "off") +
  geom_text(data = barplotData %>%
              mutate(num = ifelse(observed_expected == "observed", num , NA)), aes(label = num),
            position=position_dodge(width=0.9), vjust=-0.25, family = "Helvetica", fontface = "bold") +
  geom_text(data = barplotData %>%
              filter(observed_expected == "observed"), aes(label = paste0("p = ", pval),
                                                           y = c(29, 4, 4)),
            family = "Helvetica", color = "black", vjust = 0) +
  theme(panel.background = element_blank(),
        axis.title.x = element_blank(),
        text = element_text(family = "Helvetica"),
        legend.title = element_blank(),
        legend.position = c(0.8, 0.8),
        axis.ticks.x = element_blank(),
        axis.text = element_text(color = "black", size = 10),
        axis.line.y = element_line(linewidth = 0.25),
        axis.ticks.y = element_line(linewidth = 0.25),
        axis.text.x = element_text(face = "bold"),
        axis.ticks.length.y = unit(-0.1, "cm"))

ggsave(filename = "plots/sex_FNF_OA_barplot.pdf", width = 5, height = 5, units = "in")



# Gene scatterplots -------------------------------------------------------

# Compare to FN-f
ggplot(sex_fnf_degenes_all, aes(x = fnf_log2FoldChange, y = sex_log2FoldChange,
                            label = symbol)) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_vline(xintercept = 0, lty = 2) +
  coord_cartesian(clip = "off") +
  geom_text_repel(family = "Helvetica",
                  fontface = "bold",
                  size = 5,
                  min.segment.length = 0, seed = 3,
                  segment.color = "grey75",
                  box.padding = 0.3,
                  nudge_x = ifelse(sex_fnf_degenes_all$fnf_log2FoldChange < 0 & 
                                     !grepl("OFD1P", sex_fnf_degenes_all$symbol), -0.75, 0.75),
                  nudge_y = ifelse(sex_fnf_degenes_all$sex_log2FoldChange > 0, 
                                   ifelse(grepl("OFD1P1Y", sex_fnf_degenes_all$symbol), 0, 0.75), 
                                   -0.75)) +
  geom_point(aes(color = sex_fnf_group)) + 
  scale_x_continuous(limits = c(-6, 6), expand = c(0,0), 
                     name = '<span style = "color:#FFC154;">**FN-f** </span>log~2~(fold change)') +
  scale_y_continuous(limits = c(-10, 10), expand = c(0,0), name = "Sex log~2~(fold change)") +
  scale_color_manual(values = c("female_down" = lighten(sexColors["F"], 0.5),
                                "female_up" = darken(sexColors["F"], 0.5),
                                "male_down" = lighten(sexColors["M"], 0.5),
                                "male_up" = darken(sexColors["M"], 0.5))) +
  theme_custom_scatterplot() +
  theme(legend.position = "None",
        axis.title.x = element_markdown(size = 16),
        axis.title.y = element_markdown(size = 16),
        axis.text = element_text(color = "black", size = 12),
        axis.ticks = element_blank(),
        panel.border = element_rect(fill = NA))
ggsave(filename = "plots/sex_FNFlog2FC_scatter.pdf",
       width = 7, height = 7, units = "in")


# Compared to OA
ggplot(sex_oa_degenes_all, aes(x = oa_log2FoldChange, y = sex_log2FoldChange,
                            label = symbol)) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_vline(xintercept = 0, lty = 2) +
  coord_cartesian(clip = "off") +
  geom_text_repel(family = "Helvetica",
                  fontface = "bold",
                  size = 5,
                  min.segment.length = 0, seed = 3,
                  segment.color = NA,
                  box.padding = 0.3) +
  geom_point(aes(color = sex_oa_group)) + 
  scale_x_continuous(limits = c(-3, 3), expand = c(0,0),
                                name = '<span style = "color:#2D87BB;">**OA** </span>log~2~(fold change)') +
  scale_y_continuous(limits = c(-10, 10),
                     expand = c(0,0), name = "Sex log~2~(fold change)") +
  scale_color_manual(values = c("female_down" = lighten(sexColors["F"], 0.5),
                                "female_up" = darken(sexColors["F"], 0.5),
                                "male_down" = lighten(sexColors["M"], 0.5),
                                "male_up" = darken(sexColors["M"], 0.5))) +
  theme_custom_scatterplot() +
  theme(legend.position = "None",
        axis.title.x = element_markdown(size = 16),
        axis.title.y = element_markdown(size = 16),
        axis.text = element_text(color = "black", size = 12),
        axis.ticks = element_blank(),
        panel.border = element_rect(fill = NA))
ggsave(filename = "plots/sex_OAlog2FC_scatter.pdf",
       width = 6, height = 6, units = "in")

# OA gene violin plots ----------------------------------------------------

load("data/dds_sex_ctl.rda")
load("data/dds_sex_fnf.rda")


all_gene_sex_counts <- list()

for (geneRow in 1:nrow(sex_oa_degenes_all)){
  gene <- sex_oa_degenes_all[geneRow, ]
  
  if (gene[["condition"]] == "ctl"){
    gene_sex_counts <- get_gene_sex_Counts(gene[["gene_id"]], dds_sex_ctl)
  } else if (gene[["condition"]] == "fnf"){
    gene_sex_counts <- get_gene_sex_Counts(gene[["gene_id"]], dds_sex_fnf)
  }
  
  all_gene_sex_counts[[gene[["gene_id"]]]] <- gene_sex_counts
}

all_sex_oa_counts <- left_join(bind_rows(all_gene_sex_counts), sex_oa_degenes_all, by = "gene_id") %>%
  mutate(point_group = paste0(sex_oa_group, "_", Sex)) %>% 
  #mutate(symbol_label = paste0(symbol, "\n", str_to_title(sex_group)))
  mutate(symbol_label = ifelse(sex_group == "male", paste0('**', symbol, '**<br>',
                                                           '<span style = "color:#4788BA;">**',
                                                           str_to_title(sex_group),
                                                           '-biased**</span>'),
                               paste0('**', symbol, '**<br>',
                                      '<span style = "color:#DD8492;">**',
                                      str_to_title(sex_group),
                                      '-biased**</span>')))


all_sex_oa_counts$oa_group <- factor(all_sex_oa_counts$oa_group,
                                     levels = c("up", "down"))
all_sex_oa_counts$sex_group <- factor(all_sex_oa_counts$sex_group,
                                      levels = c("male", "female"))
all_sex_oa_counts$symbol <- factor(all_sex_oa_counts$symbol,
                                   levels = c("RARRES2", "PDPN", "TNFRSF12A", "SERPINE2"))


gene_thresholds <- all_sex_oa_counts %>%
  group_by(symbol_label) %>%
  summarise(max_sig = max(log2(count))) %>%
  left_join(all_sex_oa_counts %>% 
              dplyr::select(symbol_label, sex_log2FoldChange, oa_group) %>%
              distinct(), by = "symbol_label") %>%
  mutate(log_label = paste0("log~2~FC = ", round(sex_log2FoldChange, digits = 3)))


ggplot(all_sex_oa_counts, aes(x = log2(count), y = symbol_label, fill = Sex)) +
  geom_violin(color = NA,  alpha = 0.2) +
  scale_fill_manual(values = rep("grey25", 2)) +
  geom_point(aes(color = point_group), 
             position = position_jitterdodge(dodge.width = 0.9), size = 0.75) +
  stat_summary(position = position_dodge(width = 0.9), 
               fun = "median", geom = "crossbar", 
               width = 0.5, color = "grey25", linewidth = 0.5) +
  stat_boxplot(geom = "errorbar", position = position_dodge(width = 0.9), 
               width = 0.25, linewidth = 0.5, color = "grey25") +
  #scale_color_manual(values = c(sexColors[["F"]], sexColors[["M"]])) +
  scale_color_manual(values = c(sexColors[["F"]], "grey50", sexColors[["F"]], "grey50",
                                "grey50", sexColors[["M"]], "grey50", sexColors[["M"]])) +
  # geom_text(data = data.frame("cat" = c("Female", "Male", "Female", "Male"),
  #                             "symbol" = c("RARRES2", "PDPN", "TNFRSF12A", "SERPINE2"),
  #                             "pos" = c(25, 25, 25, 25),
  #                             "oa_group" = c("down", "down", "up", "up")),
  #           aes(x = pos, y = symbol, label = cat), 
  #           family = "Helvetica", hjust = 1,
  #           inherit.aes = FALSE) +
  ggnewscale::new_scale_fill() +
  geom_violin(color = "grey25", fill = NA, linewidth = 0.25) +
  scale_x_continuous(name = "log~2~(normalized counts)", limits = c(0, 25), 
                     expand = c(0, 0), breaks = seq(5, 25, 5)) +
  coord_cartesian(clip = "off") +
  facet_wrap(~oa_group, ncol = 1, scales = "free_y", strip.position = "left",
             labeller = as_labeller(c("up" = "**Upregulated** in OA", "down" = "**Downregulated** in OA"))) +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        text = element_text(family = "Helvetica"),
        plot.margin = margin(r = 10, t = 10),
        strip.background = element_blank(),
        strip.placement = "outside",
        strip.text.y.left =  element_markdown(family = "Helvetica", size = 12),
        legend.position = "None",
        axis.ticks = element_blank(),
        panel.spacing = unit(0, "mm"),
        axis.title.y = element_blank(),
        axis.title.x = element_markdown(size = 12),
        axis.line = element_line(),
        axis.text.x = element_text(color = "black", size = 12),
        axis.text.y = element_markdown(color = "black", size = 12, family = "Helvetica", halign = 0,
                                       hjust = 0))

ggsave(filename = "plots/sex_OA_count_violins.pdf", width = 5, height = 6, units = "in")



## Different orientation 

ggplot(all_sex_oa_counts, aes(x = symbol_label, y = log2(count), fill = Sex)) +
  geom_violin(color = NA,  alpha = 0.2) +
  scale_fill_manual(values = rep("grey25", 2)) +
  geom_point(aes(color = point_group), 
             position = position_jitterdodge(dodge.width = 0.9), size = 0.75) +
  geom_richtext(data = gene_thresholds, aes(y = max_sig + 0.5, label = log_label, fill = NA,
                                            group = oa_group),
                family = "Helvetica", fill = NA, label.size = 0, hjust = 0.5, vjust = 0,
                label.color = NA) +
  stat_summary(position = position_dodge(width = 0.9), 
               fun = "median", geom = "crossbar", 
               width = 0.5, color = "grey25", linewidth = 0.5) +
  stat_boxplot(geom = "errorbar", position = position_dodge(width = 0.9), 
               width = 0.25, linewidth = 0.5, color = "grey25") +
  # scale_color_manual(values = c(sexColors[["F"]], "grey50", sexColors[["F"]], "grey50",
  #                               "grey50", sexColors[["M"]], "grey50", sexColors[["M"]])) +
  scale_color_manual(values = c(sexColors[["F"]], sexColors[["M"]], 
                                sexColors[["F"]], sexColors[["M"]],
                                sexColors[["F"]], sexColors[["M"]], 
                                sexColors[["F"]], sexColors[["M"]])) + 
  ggnewscale::new_scale_fill() +
  geom_violin(color = "grey25", fill = NA, linewidth = 0.25) +
  scale_y_continuous(name = "log~2~(normalized counts)", limits = c(0, 25), 
                     expand = c(0, 0), breaks = seq(5, 25, 5)) +
  coord_cartesian(clip = "off") +
  facet_wrap(~oa_group, ncol = 2, scales = "free_x", strip.position = "top",
             labeller = as_labeller(c("up" = "**Upregulated** in OA", "down" = "**Downregulated** in OA"))) +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        text = element_text(family = "Helvetica"),
        plot.margin = margin(r = 10, t = 10),
        strip.background = element_blank(),
        strip.placement = "outside",
        strip.text.x.top = element_markdown(family = "Helvetica", size = 12),
        legend.position = "None",
        panel.spacing = unit(0, "mm"),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(color = "black"),
        axis.ticks.length.y = unit(-0.1, "cm"),
        axis.title.x = element_blank(),
        axis.title.y = element_markdown(size = 12),
        axis.line = element_line(),
        axis.text.y = element_text(color = "black", size = 12),
        axis.text.x = element_markdown(color = "black", size = 12, family = "Helvetica", halign = 0.5,
                                       hjust = 0.5))

ggsave(filename = "plots/sex_OA_count_violins_horizontal.pdf", width = 8, height = 5, units = "in")

# OA gene counts in PBS and FNF -------------------------------------------




load("data/dds_sex_ctl.rda")
load("data/dds_sex_fnf.rda")

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


ctl_sex_OAgenes <- read_csv("data/ctl_sexDE_pval01.csv", col_select = c("gene_id", "symbol", "padj")) %>%
  filter(gene_id %in% sex_oa_degenes_all$gene_id) %>%
  dplyr::rename(PBS = padj)
fnf_sex_OAgenes <- read_csv("data/fnf_sexDE_pval01.csv", col_select = c("gene_id", "symbol", "padj")) %>%
  filter(gene_id %in% sex_oa_degenes_all$gene_id) %>%
  dplyr::rename(`FN-f` = padj)

all_sig_sex_OAgenes <- full_join(ctl_sex_OAgenes, fnf_sex_OAgenes) %>%
  pivot_longer(cols = c("PBS", "FN-f"), names_to = "condition", values_to = "padj") %>%
  left_join(gene_maxCounts) %>%
  filter(!is.na(padj)) %>%
  mutate(signif = "*")
all_sig_sex_OAgenes <- merge(all_sig_sex_OAgenes, data.frame("Sex" = c("M", "F")))




up_boxplots <- ggplot(all_sex_oa_counts %>% filter(oa_group == "up"), aes(x = Sex, y = log2(count), fill = Sex)) +
  geom_boxplot() +
  ggh4x::facet_grid2(condition ~ symbol, switch = "both", scales = "free",
                     independent = "all") +
  geom_line(data = all_sig_sex_OAgenes %>% filter(symbol %in% c("SERPINE2", "TNFRSF12A")), aes(x = Sex, y = maxCount + 0.5, group = symbol), inherit.aes = FALSE,
            linewidth = 0.5) +
  geom_text(data = all_sig_sex_OAgenes %>% filter(symbol %in% c("SERPINE2", "TNFRSF12A")), aes(x = 1.5, y = Inf, label = signif), inherit.aes = FALSE,
            fontface = "bold") +
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


down_boxplots <- ggplot(all_sex_oa_counts %>% filter(oa_group == "down"), aes(x = Sex, y = log2(count), fill = Sex)) +
  geom_boxplot() +
  ggh4x::facet_grid2(condition ~ symbol, switch = "both", scales = "free",
                     independent = "all") +
  geom_line(data = all_sig_sex_OAgenes %>% filter(symbol %in% c("PDPN", "RARRES2")), aes(x = Sex, y = maxCount + 0.5, group = symbol), inherit.aes = FALSE,
            linewidth = 0.5) +
  geom_text(data = all_sig_sex_OAgenes %>% filter(symbol %in% c("PDPN", "RARRES2")), aes(x = 1.5, y = Inf, label = signif), inherit.aes = FALSE,
            fontface = "bold") +
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


up_boxplots + down_boxplots


ggsave(filename = "plots/oa_sex_genes_boxplots.pdf", width = 9, height = 5, units= "in")

