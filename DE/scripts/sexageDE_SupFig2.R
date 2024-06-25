library(tidyverse)
source("../plotting_utils.R")

# Plotting function -------------------------------------------------------

plotSexResponseGene <- function(geneRow, dds){
  
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
    scale_color_manual(values = c(sexColors[["F"]], sexColors[["M"]])) +
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
  ggsave(plot, filename = paste0("plots/sexDE_Fig2_supp/sexspecifictreatment_",
                                 genesymbol, ".pdf"),
         width = 6, height = 6, units = "in" )
  return(plot)
  
} 

# Load and plot genes -----------------------------------------------------
load("data/sex_de/dds_sex_treatment_response.rda")
sexDE_treatmenteffect_pval01 <- 
  read_csv("data/sex_de/sexDE_treatmenteffect_pval01.csv")

apply(sexDE_treatmenteffect_pval01, 
      1, plotSexResponseGene, dds = dds_sex_treatment_response)

# CHROMOSOME BARPLOTS FOR SEX DE GENES -----------------------------------------

ctl_sig_genes <- read_csv("data/sex_de/ctl_sexDE_pval01.csv") |> 
  mutate(condition = "PBS") |> 
  mutate(sex = ifelse(log2FoldChange < 0, "Female", "Male"))
fnf_sig_genes <- read_csv("data/sex_de/fnf_sexDE_pval01.csv") |> 
  mutate(condition = "FN-f") |> 
  mutate(sex = ifelse(log2FoldChange < 0, "Female", "Male"))

all_sig_sex_genes <- bind_rows(ctl_sig_genes,
                               fnf_sig_genes)
all_sig_sex_genes$seqnames <- factor(all_sig_sex_genes$seqnames,
                                     levels = c(as.character(1:22), "X", "Y"))
all_sig_sex_genes$condition <- factor(all_sig_sex_genes$condition,
                                      levels = c("PBS", "FN-f"))

sex_fnfpbs_numGenes_barplot <- ggplot(all_sig_sex_genes, aes(x = seqnames, fill = sex)) +
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
        panel.background = element_rect(fill = "transparent", color = "transparent"),
        plot.background = element_rect(fill = "transparent", color = "transparent"),
        panel.spacing = unit(0, "cm"),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(linewidth = 0.25),
        axis.ticks.length.y = unit(-0.1, "cm"),
        axis.title.y = element_markdown(family = "Helvetica", size = 7),
        axis.title.x = element_text(size = 8),
        axis.line = element_line(linewidth = 0.25),
        axis.text.x = element_text(color = "black", size = 6,
                                   face = c(rep("plain", 22), rep("bold", 2))),
        axis.text.y = element_text(color = "black", size = 6),
        strip.background = element_blank(),
        strip.text = element_text(size = 8, vjust = 0, face = "bold"))

save(sex_fnfpbs_numGenes_barplot, file = "plots/sexageDE_Fig2_supp/sex_fnfpbs_numGenes_barplot.rda")


# VENN DIAGRAM OF PBS AND FNF SEX GENES -----------------------------------
ctl_sig_genes <- read_csv("data/sex_de/ctl_sexDE_pval01.csv")
fnf_sig_genes <- read_csv("data/sex_de/fnf_sexDE_pval01.csv")

pbs_fnf_sex_genes <- tibble(values = unique(c(ctl_sig_genes$gene_id, fnf_sig_genes$gene_id))) |> 
  mutate(PBS = values %in% ctl_sig_genes$gene_id,
         FNF = values %in% fnf_sig_genes$gene_id)

sex_pbs_fnf_venndiagram <- ggplot(pbs_fnf_sex_genes, aes(A = PBS, B = FNF)) +
  geom_venn(set_names = c("PBS", "FN-f"), 
            fill_color = c(log2fcColors[["-"]], log2fcColors[["+"]]), 
            stroke_color = NA, auto_scale = TRUE, show_percentage = FALSE,
            text_size = 3, set_name_size = 3) +
  coord_fixed()  +
  theme(panel.background = element_rect(fill = "transparent", color = "transparent"),
        plot.background = element_rect(fill = "transparent", color = "transparent"),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank())
sex_pbs_fnf_venndiagram <- venn_font(sex_pbs_fnf_venndiagram, font = "Helvetica")
save(sex_pbs_fnf_venndiagram, file = "plots/sexageDE_Fig2_supp/sex_pbs_fnf_venndiagram.rda")

#### Checking directions of effect for each set

# Overlap
overlap_sex_genes <- ctl_sig_genes[which(ctl_sig_genes$gene_id %in% fnf_sig_genes$gene_id), "gene_id"]
ctl_overlap <- ctl_sig_genes |> 
  filter(gene_id %in% overlap_sex_genes$gene_id) |> 
  dplyr::select(gene_id, symbol, log2FoldChange) |> 
  dplyr::rename(ctl_l2fc = log2FoldChange)

fnf_overlap <- fnf_sig_genes |> 
  filter(gene_id %in% overlap_sex_genes$gene_id) |> 
  dplyr::select(gene_id, symbol, log2FoldChange) |> 
  dplyr::rename(fnf_l2fc = log2FoldChange)

overlap_directions <- left_join(ctl_overlap,
                                fnf_overlap,
                                by = c("gene_id", "symbol")) |> 
  mutate(concordant_dir = ifelse(sign(ctl_l2fc) == sign(fnf_l2fc), TRUE, FALSE))

overlap_percent_concordant <- length(which(overlap_directions$concordant_dir))/nrow(overlap_directions)

# ctl only, looking up in all fnf
ctl_unique <- ctl_sig_genes[which(!ctl_sig_genes$gene_id %in% fnf_sig_genes$gene_id), "gene_id"]

ctl_ctl_unique <- ctl_sig_genes |> 
  filter(gene_id %in% ctl_unique$gene_id) |> 
  dplyr::select(gene_id, symbol, log2FoldChange) |> 
  dplyr::rename(ctl_l2fc = log2FoldChange)

fnf_ctl_unique <- read_csv("data/sex_de/fnf_sex_shrink.csv",
                           col_select = c("gene_id", "symbol", "log2FoldChange")) |> 
  filter(gene_id %in% ctl_unique$gene_id) |> 
  dplyr::rename(fnf_l2fc = log2FoldChange)

ctl_unique_directions <- left_join(ctl_ctl_unique,
                                   fnf_ctl_unique,
                                   by = c("gene_id", "symbol")) |> 
  mutate(concordant_dir = ifelse(sign(ctl_l2fc) == sign(fnf_l2fc), TRUE, FALSE))

ctl_unique_percent_concordant <- length(which(ctl_unique_directions$concordant_dir))/nrow(ctl_unique_directions)

# fnf only, looking up in all ctl
fnf_unique <- fnf_sig_genes[which(!fnf_sig_genes$gene_id %in% ctl_sig_genes$gene_id), "gene_id"]

fnf_fnf_unique <- fnf_sig_genes |> 
  filter(gene_id %in% fnf_unique$gene_id) |> 
  dplyr::select(gene_id, symbol, log2FoldChange) |> 
  dplyr::rename(fnf_l2fc = log2FoldChange)

ctl_fnf_unique <- read_csv("data/sex_de/ctl_sex_shrink.csv",
                           col_select = c("gene_id", "symbol", "log2FoldChange")) |> 
  filter(gene_id %in% fnf_unique$gene_id) |> 
  dplyr::rename(ctl_l2fc = log2FoldChange)

fnf_unique_directions <- left_join(fnf_fnf_unique,
                                   ctl_fnf_unique,
                                   by = c("gene_id", "symbol")) |> 
  mutate(concordant_dir = ifelse(sign(ctl_l2fc) == sign(fnf_l2fc), TRUE, FALSE))

fnf_unique_percent_concordant <- length(which(fnf_unique_directions$concordant_dir))/nrow(fnf_unique_directions)

# VENN DIAGRAM OF PBS AND FNF AGE GENES -----------------------------------

ctl_cluster_pval05 <- read_csv("data/age_de/ctl_age_pval05clusters.csv") %>%
  distinct(gene_id, .keep_all = TRUE) %>%
  dplyr::select(gene_id, cluster)

fnf_cluster_pval05 <- read_csv("data/age_de/fnf_age_pval05clusters.csv") %>%
  distinct(gene_id, .keep_all = TRUE) %>%
  dplyr::select(gene_id, cluster)

pbs_fnf_age_genes <- tibble(values = unique(c(ctl_cluster_pval05$gene_id, fnf_cluster_pval05$gene_id))) %>%
  mutate(PBS = values %in% ctl_cluster_pval05$gene_id,
         FNF = values %in% fnf_cluster_pval05$gene_id)

age_pbs_fnf_venndiagram <- ggplot(pbs_fnf_age_genes, aes(A = PBS, B = FNF)) +
  geom_venn(set_names = c("PBS", "FN-f"), 
            fill_color = c(log2fcColors[["-"]], log2fcColors[["+"]]), 
            stroke_color = NA, auto_scale = TRUE, show_percentage = FALSE,
            text_size = 4, set_name_size = 4) +
  coord_fixed()  +
  theme(panel.background = element_rect(fill = "transparent", color = "transparent"),
        plot.background = element_rect(fill = "transparent", color = "transparent"),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank())
age_pbs_fnf_venndiagram <- venn_font(age_pbs_fnf_venndiagram, font = "Helvetica")

save(age_pbs_fnf_venndiagram, file = "plots/sexageDE_Fig2_supp/age_pbs_fnf_venndiagram.rda")




# AGE OVERLAP WITH AGE-RELATED GENES FROM GTEX ----------------------------

# Get union of significant age genes
ctl_cluster_pval05 <- read_csv("data/age_de/ctl_age_pval05clusters.csv") |> 
  distinct(gene_id, .keep_all = TRUE) |> 
  dplyr::select(gene_id, symbol, cluster) 

fnf_cluster_pval05 <- read_csv("data/age_de/fnf_age_pval05clusters.csv") |> 
  distinct(gene_id, .keep_all = TRUE) |> 
  dplyr::select(gene_id, symbol, cluster)

union_sig_genes <- union(ctl_cluster_pval05$gene_id, 
                         fnf_cluster_pval05$gene_id)

age_genes <- bind_rows(ctl_cluster_pval05 |> 
                         filter(gene_id %in% union_sig_genes),
                       fnf_cluster_pval05 |> 
                         filter(gene_id %in% union_sig_genes)) |> 
  distinct(gene_id, .keep_all = TRUE) |> 
  mutate(tissue = "Chondrocytes")

gtex_age_de <- read_csv("data/Yang_age/41598_2015_BFsrep15145_MOESM2_ESM.csv") |> 
  dplyr::rename(gene_id = `Gene ID`) |> 
  dplyr::select(-ends_with("p-value")) |> 
  # convert IDs to ones compatible with ours
  mutate(gene_id = gsub("\\..*", "", gene_id)) |> 
  # pivot longer with different tissues
  pivot_longer(cols = c(-gene_id, -`Gene Symbol`,- `Entrez ID`),
               names_to = "name",
               values_to = "value") |> 
  separate_wider_delim(cols = "name",delim = "_", names = c("tissue", "stat")) |> 
  pivot_wider(names_from = "stat", values_from = "value") |> 
  filter(FDR < 0.05) |> 
  mutate(cluster = ifelse(AgeCoef > 0, "+", "-")) |> 
  # filter for our genes
  filter(gene_id %in% age_genes$gene_id) |> 
  # Omit Skin and Thyroid like paper
  filter(!tissue %in% c("Skin", "Thyroid"))
# Omit skin and thyroid like paper?

# Join with our data
all_data_age_genes <- bind_rows(age_genes  |> 
                                  dplyr::select(gene_id, tissue, cluster),
                                gtex_age_de |>  
                                  dplyr::select(gene_id, tissue, cluster)) |> 
  complete(gene_id, tissue) |> 
  left_join(age_genes |> dplyr::select(gene_id, symbol), by = "gene_id") 

# Get order of genes by number of overlaps between tissues
age_geneOverlaps <- all_data_age_genes |> 
  group_by(symbol) |> 
  filter(!is.na(cluster)) |> 
  summarize(nOverlap = dplyr::n()) |> 
  arrange(nOverlap)

# Join gene ordering with merged gtex/pbs/fnf sb gene data
all_data_age_genes <- left_join(all_data_age_genes, age_geneOverlaps, by = "symbol") |> 
  group_by(nOverlap) |> 
  # Make all nOverlap groups consistent with age cluster ordering
  arrange(cluster, .by_group = TRUE) |> 
  ungroup()

# Set factors so chondrocytes are in the first column and gene symbols are ordered
# by number of overlaps
all_data_age_genes$tissue = factor(all_data_age_genes$tissue, levels = c("Chondrocytes", 
                                                                         unique(sort(gtex_age_de$tissue))))
all_data_age_genes$symbol <- factor(all_data_age_genes$symbol, levels = unique(all_data_age_genes$symbol))

all_data_age_genes <- all_data_age_genes |> 
  mutate(symbol_order_temp = factor(symbol, levels = rev(levels(symbol)))) |> 
  # Assign gene symbol and tissues numbers based on factor ordering
  mutate(tissue_no = as.numeric(tissue),
         gene_no = as.numeric(symbol_order_temp)) |> 
  dplyr::select(-symbol_order_temp)

gtex_age <- ggplot(all_data_age_genes, aes(x = tissue, y = symbol)) +
  scale_x_discrete(position = "top") +
  coord_cartesian(clip = "off") +
  geom_tile(aes(fill = cluster), color = "grey20", linewidth = 0.005) +
  geom_rect(data = tibble(xmin = 0.5, xmax = 1.5, ymin = 0.5,
                          ymax = length(levels(all_data_age_genes$symbol))+0.5),
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = NA, color = "black", inherit.aes = FALSE, linewidth = 0.25) +
  geom_rect(data = tibble(xmin = 0.5, xmax = length(levels(all_data_age_genes$tissue)) + 0.5,
                          ymin = 0.5, 
                          ymax = (max(all_data_age_genes$gene_no)+1 - all_data_age_genes |> 
                                    filter(nOverlap == 1) |> 
                                    pull(gene_no) |> 
                                    min()) + 0.5),
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = NA, color = "black", inherit.aes = FALSE, linewidth = 0.25) +
  scale_fill_manual(values = c("-" = ageClusterColors[["-"]], 
                               "+" = ageClusterColors[["+"]]), na.value = NA) +
  scale_y_discrete(expand = c(0, 0)) +
  theme(panel.background = element_rect(fill = "transparent", color = "transparent"),
        plot.background = element_rect(fill = "transparent", color = "transparent"),
        axis.text.x = element_text(family = "Helvetica", color= "black", size = 5, 
                                   angle = 20, hjust = 0.5),
        axis.text.y = element_blank(),
        axis.title = element_blank(),
        text = element_text(family = "Helvetica"),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(linewidth = 0.2),
        axis.ticks.length = unit(1.5, "mm"),
        legend.position = "none")

save(gtex_age, file = "plots/sexageDE_Fig2_supp/gtex_age.rda")

# Assemble supp figure with pg --------------------------------------------

pdf("plots/sexageDE_Fig2_supp/SupFig2.pdf",
    width = 7, height = 4.1)
pageCreate(width = 7, height = 4.1, showGuides = FALSE)
# A - Sex chromosome barplots
plotText("A", x = 0.1, y = 0.1, just = c("left", "top"), fontfamily = "Helvetica",
         fontsize = 14, fontface = "bold")
plotGG(sex_fnfpbs_numGenes_barplot, x = 0.1, y = 0.1, width = 4, height = 2.25)

# B - Sex Venn Diagram
plotText("B", x = 0.1, y = 2.4, just = c("left", "top"), fontfamily = "Helvetica",
         fontsize = 14, fontface = "bold")
plotGG(sex_pbs_fnf_venndiagram, x = -0.15, y = 2.1, width = 2.5, height = 2.5)
plotText(label = "Sex-specific genes",
         fontfamily = "Helvetica", fontsize = 8, x = 1, y = 3.95)

# C - Age Venn Diagram
plotText("C", x = 2, y = 2.4, just = c("left", "top"), fontfamily = "Helvetica",
         fontsize = 14, fontface = "bold")
plotGG(age_pbs_fnf_venndiagram, x = 1.7, y = 2.1, width = 2.5, height = 2.5)
plotText(label = "Age-related genes",
         fontfamily = "Helvetica", fontsize = 8, x = 3, y = 3.95)

# D - age GTEx
plotText("D", x = 4.05, y = 0.1, just = c("left", "top"), fontfamily = "Helvetica",
         fontsize = 14, fontface = "bold")

plotGG(gtex_age, x = 4.2, y = -0.1, width = 2.75, height = 4.2)

dev.off()

# Make supplemental table of sex de genes with FNF and OA info ------------


# FN-f DE
fnf_de_genes <- read_csv("data/condition_de/sig_deGenes_pval05_l2fc1.csv",
                         col_select = c("gene_id", "log2FoldChange", "padj")) |>
  mutate(`FNF response` = case_when(padj < 0.01 &
                                      abs(log2FoldChange) > 2 &
                                      log2FoldChange < 0 ~ "---",
                                    padj < 0.01 &
                                      abs(log2FoldChange) > 2 &
                                      log2FoldChange > 0 ~ "+++", 
                                    padj < 0.05 &
                                      abs(log2FoldChange) > 1 &
                                      log2FoldChange < 0 ~ "-",
                                    padj < 0.05 & 
                                      abs(log2FoldChange) > 1 &
                                      log2FoldChange > 0 ~ "+")) |> 
  dplyr::select(-padj, -log2FoldChange)

# OA DE 
sex_oa_degenes <- read_csv("data/sex_de/sex_oa_degenes_all.csv",
                           col_select = c("gene_id", "oa_group", "oa_study")) |> 
  dplyr::rename(`Change with OA` = oa_group,
                `OA study` = oa_study)

# Union of PBS and FNF results
ctl_sex_degenes <- read_csv("data/sex_de/ctl_sexDE_pval01.csv", 
                            col_select = c("symbol", "gene_id", "padj", "log2FoldChange")) |> 
  mutate(`M/F` = ifelse(log2FoldChange < 0, "F", "M")) |> 
  dplyr::rename(sex_padj_PBS = padj,
                sex_log2FoldChange_PBS = log2FoldChange)
fnf_sex_degenes <- read_csv("data/sex_de/fnf_sexDE_pval01.csv",
                            col_select = c("symbol", "gene_id", "padj", "log2FoldChange")) |> 
  mutate(`M/F` = ifelse(log2FoldChange < 0, "F", "M")) |> 
  dplyr::rename(sex_padj_FNF = padj,
                sex_log2FoldChange_FNF = log2FoldChange)


union_sex_degenes <- full_join(ctl_sex_degenes, 
                               fnf_sex_degenes, by = c("symbol", "gene_id", "M/F")) |> 
  relocate("M/F", .after = sex_log2FoldChange_FNF) |> 
  left_join(fnf_de_genes, by = "gene_id") |> 
  left_join(sex_oa_degenes, by = "gene_id") |> 
  arrange(symbol)

write_csv(union_sex_degenes, file = "tables/SupTable3.csv")

# Write to google drive
ss <- gs4_create(name = "SupTable3")
write_sheet(union_sex_degenes,
            ss, sheet = "Sheet1")


drive_mv(file = "SupTable3", path = as_dribble("CQTL paper/Figures and Tables"))

# Make supplemental table of age de genes with FNF and OA info ---------------------------------------

# FN-f DE
fnf_de_genes <- read_csv("data/condition_de/sig_deGenes_pval05_l2fc1.csv",
                         col_select = c("gene_id", "log2FoldChange", "padj")) |>
  mutate(`FNF response` = case_when(padj < 0.01 &
                                      abs(log2FoldChange) > 2 &
                                      log2FoldChange < 0 ~ "---",
                                    padj < 0.01 &
                                      abs(log2FoldChange) > 2 &
                                      log2FoldChange > 0 ~ "+++", 
                                    padj < 0.05 &
                                      abs(log2FoldChange) > 1 &
                                      log2FoldChange < 0 ~ "-",
                                    padj < 0.05 & 
                                      abs(log2FoldChange) > 1 &
                                      log2FoldChange > 0 ~ "+")) |> 
  dplyr::select(-padj, -log2FoldChange)

## OA
age_oa_degenes <- read_csv("data/age_de/age_oa_degenes_all.csv",
                           col_select = c("gene_id", "oa_group", "oa_study")) |> 
  dplyr::rename(`Change with OA` = oa_group,
                `OA study` = oa_study)

pbs_age <- read_csv("data/age_de/ctl_age_pval05clusters.csv", col_select = c("symbol", "gene_id", "padj", "cluster")) |> 
  distinct() |> 
  dplyr::rename(`Change with age` = cluster,
                age_padj_PBS = padj) |> 
  mutate(`Change with age` = ifelse(`Change with age` == "-", "down", "up"))

fnf_age <- read_csv("data/age_de/fnf_age_pval05clusters.csv", col_select = c("symbol", "gene_id", "padj", "cluster")) |> 
  distinct() |> 
  dplyr::rename(`Change with age` = cluster,
                age_padj_FNF = padj) |> 
  mutate(`Change with age` = ifelse(`Change with age` == "-", "down", "up"))

pbs_fnf_age_table <- full_join(pbs_age, fnf_age, by = c("symbol", "gene_id", "Change with age")) |> 
  relocate("Change with age", .after = age_padj_FNF) |> 
  left_join(fnf_de_genes, by = "gene_id") |> 
  left_join(age_oa_degenes, by = "gene_id") |> 
  group_by(symbol) |> 
  mutate(`Change with age` = paste(`Change with age`, collapse = "/"),
         age_padj_PBS = paste(na.omit(age_padj_PBS), collapse = "/"),
         age_padj_FNF = paste(na.omit(age_padj_FNF), collapse = "/")) |> 
  ungroup() |> 
  distinct() |> 
  mutate(age_padj_PBS = as.numeric(age_padj_PBS),
         age_padj_FNF = as.numeric(age_padj_FNF)) |> 
  arrange(symbol)

write_csv(pbs_fnf_age_table, file = "tables/SupTable5.csv") 

# Write to google drive
ss <- gs4_create(name = "SupTable5")
write_sheet(pbs_fnf_age_table,
            ss, sheet = "Sheet1")
drive_mv(file = "SupTable5", path = as_dribble("CQTL paper/Figures and Tables"))

# Assemble Age GO/KEGG Supp table ---------------------------------------------

upage_go <- read_csv("tables/SupTable6A.csv")
downage_go <- read_csv("tables/SupTable6B.csv")


ss <- gs4_create(name = "SupTable6")
write_sheet(upage_go,
            ss, sheet = "GO Terms - Up with Age")
write_sheet(downage_go,
            ss, sheet = "GO Terms - Down with Age")

drive_mv(file = "SupTable6", path = as_dribble("CQTL paper/Figures and Tables"))
