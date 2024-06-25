library(DESeq2)
library(tidyverse)
library(DescTools)
library(ggh4x)
library(scales)
library(patchwork)
library(googledrive)
library(googlesheets4)
source("../plotting_utils.R")


# Functions ---------------------------------------------------------------

cor_amongDonors <- function(geneCountMatrix, condition){
  
  # Get donors based on column names
  donors <- unique(unlist(lapply(str_split(colnames(geneCountMatrix), "_"), 
                          `[[`, 2)))
  amongCorrelations <- c()
  for (donor in donors){
    
    # Pull out the paired replicate donor gene counts
    geneCountDonor <- geneCountMatrix |> 
      dplyr::select(contains(donor))
    
    # Calculate Pearson's correlation coefficient
    pearsons_coef <- cor(geneCountDonor[,1], geneCountDonor[,2], 
                         method = "pearson")
    amongCorrelations[donor] <- pearsons_coef
  }
  
  # Make dataframe
  amongDonors_df <- data.frame("coef" = amongCorrelations,
                               "Donor" = names(amongCorrelations)) |> 
    remove_rownames() |> 
    mutate(group = "within",
           Condition = condition,
           fisherz = FisherZ(coef))
  
  return(amongDonors_df)
}

cor_acrossDonors <- function(repCountMatrix, standCountMatrix, condition){
  
  # Get replicate donors based on column names
  repdonors <- unique(unlist(lapply(str_split(colnames(repCountMatrix), "_"), 
                                 `[[`, 2)))
  acrossCorrelations <- list()
  for (donor in repdonors){
    
    # Pull out 2nd donor replicate data 
    repCountDonor <- repCountMatrix |> 
      dplyr::select(contains(donor)) |> 
      dplyr::select(ends_with("_2")) |> 
      rownames_to_column(var = "gene_id")
    
    # Filter original donor data from standard count matrix
    otherDonors <- standCountMatrix |> 
      dplyr::select(!contains(donor)) |> 
      rownames_to_column(var = "gene_id")
    
    # Join to get same genes
    fullMatrix <- inner_join(repCountDonor, otherDonors, by = "gene_id") |> 
      column_to_rownames(var = "gene_id")
    
    # Calculate Pearson's correlation coefficients
    pearsons_coefs <- cor(fullMatrix[,1], fullMatrix[,-1], 
                         method = "pearson")
    
    acrossCorrelations[[donor]] <- data.frame("coef" = pearsons_coefs[1,],
                                              "Donor" = donor)
  }
  
  # Make dataframe of all results
  acrossDonors_df <- bind_rows(acrossCorrelations) |> 
    remove_rownames() |> 
    mutate(group = "across",
           Condition = condition,
           fisherz = FisherZ(coef))
  
  return(acrossDonors_df)
  
}

# Library read depth barplot ---------------------------------------------------

mqc_cutadapt_filtered_reads <- read_delim("data/qc/mqc_cutadapt_filtered_reads_plot_1.txt") |> 
  separate_wider_delim(cols = "Sample", delim = "_", names = c(NA, "Donor", NA, 
                                                               "Condition", NA, 
                                                               NA, "Read"),
                       cols_remove = FALSE) |> 
  filter(Read == "R1") |> 
  mutate(Sample = gsub("_R1", "", Sample)) |> 
  # Remove geno-contaminated donor and donor with poor STAR alignment
  filter(Donor != "AM7352" & Donor != "AM7244")

mean_depth <- mean(mqc_cutadapt_filtered_reads$`Reads passing filters`)

library_read_depth_barplot <- ggplot(mqc_cutadapt_filtered_reads, aes(x = Sample, 
                                        y = `Reads passing filters`)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = mean_depth, lty = 2, linewidth = 0.25) +
  scale_x_discrete(name = "Library", expand = expansion(mult = c(0.005, 0))) +
  scale_y_continuous(expand = c(0, 0),
                     labels = label_number(scale = 1e-6, big.mark = ""),
                     breaks = seq(0, 150e6, 25e6),
                     name = "Read depth (millions of reads)") +
  theme_custom_general() + 
  theme(text = element_text(family = "Helvetica"),
        axis.text.x = element_blank(),
        axis.title = element_text(family = "Helvetica", size = 8))

save(library_read_depth_barplot, 
     file = "plots/conditionDE_Fig1_supp/library_read_depth_barplot.rda")

# Violin plots of QCs ------------------------------------------------------

# Library read depth
mqc_cutadapt_filtered_reads <- read_delim("data/qc/mqc_cutadapt_filtered_reads_plot_1.txt") |> 
  separate_wider_delim(cols = "Sample", delim = "_", names = c(NA, "Donor", NA, 
                                                               "Condition", NA, 
                                                               NA, "Read"),
                       cols_remove = FALSE) |> 
  filter(Read == "R1") |> 
  mutate(Sample = gsub("_R1", "", Sample)) |> 
  # Remove geno-contaminated donor and donor with poor STAR alignment
  filter(Donor != "AM7352" & Donor != "AM7244") |> 
  mutate(group = "Library read depths")


read_depth_violin <- ggplot(mqc_cutadapt_filtered_reads, aes(x = group, y = `Reads passing filters`)) +
  geom_violin(color = NA, fill = "#78A1Cd", alpha = 0.6) +
  geom_violin(fill = NA) +
  stat_boxplot(geom = "errorbar", width = 0.1, linewidth = 0.75) +
  stat_summary(fun = "mean", geom = "crossbar", width = 0.1, linewidth = 0.5) +
  scale_y_continuous(expand = c(0, 0),
                     labels = label_number(scale = 1e-6, big.mark = ""),
                     name = "Read depth (millions of reads)",
                     limits = c(60e6, 150e6)) +
  theme_custom_scatterplot() +
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(family = "Helvetica", color = "black",
                                   size = 8, vjust = 0),
        axis.title.y = element_text(family = "Helvetica", color = "black",
                                   size = 8))


# Library RIN values
library_RINs <- read_csv("data/samplesheet.csv", col_select = c("Sample", "Donor", "RIN")) |> 
  filter(Donor != "AM7352") |> 
  distinct() |> 
  mutate(group = "Library RIN\nvalues")

rin_violin <- ggplot(library_RINs, aes(x = group, y = RIN)) +
  geom_violin(color = NA, fill = "#A9D48A", alpha = 0.6) +
  geom_violin(fill = NA) +
  stat_boxplot(geom = "errorbar", width = 0.1, linewidth = 0.75) +
  stat_summary(fun = "mean", geom = "crossbar", width = 0.1, linewidth = 0.5) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(8.5, 10.1)) +
  theme_custom_scatterplot() +
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(family = "Helvetica", color = "black",
                                   size = 8, vjust = 0),
        axis.title.y = element_text(family = "Helvetica", color = "black",
                                   size = 8))

# Salmon mapping
salmon_mapping <- read_delim("data/qc/multiqc_salmon.txt", col_select = c("Sample",
                                                                          "percent_mapped")) |> 
  separate_wider_delim(cols = "Sample", delim = "_", names = c(NA, "Donor", NA, 
                                                               "Condition", NA, 
                                                               NA),
                       cols_remove = FALSE) |> 
  # Remove geno-contaminated donor and donor with poor STAR alignment
  filter(Donor != "AM7352" & Donor != "AM7244") |> 
  mutate(group = "Library mapping\nwith Salmon")

salmon_mapping_violin <- ggplot(salmon_mapping, aes(x = group, y = percent_mapped)) +
  geom_violin(color = NA, fill = "#FBBE67", alpha = 0.6) +
  geom_violin(fill = NA) +
  stat_boxplot(geom = "errorbar", width = 0.1, linewidth = 0.75) +
  stat_summary(fun = "mean", geom = "crossbar", width = 0.1, linewidth = 0.5) +
  scale_y_continuous(expand = c(0, 0), name = "Percent mapped (%)",
                     limits = c(85, 100)) +
  theme_custom_scatterplot() +
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(family = "Helvetica", color = "black",
                                   size = 8, vjust = 0),
        axis.title.y = element_text(family = "Helvetica", color = "black",
                                   size = 8))

qc_violins <- read_depth_violin + 
  rin_violin + salmon_mapping_violin + 
  plot_annotation(theme = theme(panel.background = 
                                  element_rect(fill = "transparent", 
                                               color = "transparent"),
                                plot.background = 
                                  element_rect(fill = "transparent", 
                                               color = "transparent")))

save(qc_violins, file = "plots/conditionDE_Fig1_supp/qc_violins.rda")

# PCA ---------------------------------------------------------------------

load("data/condition_de/differential_expression_dds.rda")
donorSamplesheet <- read_csv("data/donorSamplesheet.csv")

# Normalized gene counts
normCounts <- t(assay(vst(dds)))

# Calculate principal components and percent variance explained
pcs_normCounts <- prcomp(normCounts)
percentVar <- pcs_normCounts$sdev^2/sum(pcs_normCounts$sdev^2)

# Create dataframe for plotting
pc_df <- data.frame("PC1" = pcs_normCounts$x[,1],
                    "PC2" = pcs_normCounts$x[,2]) %>%
  rownames_to_column(var = "Sample")

# Compile donor sex and sample treatment for coloring covariates
covariates <- as.data.frame(colData(dds)[,c("Condition", "Donor")]) %>%
  rownames_to_column(var = "Sample") %>%
  left_join(donorSamplesheet[,c("Donor", "Sex")], by = "Donor")

pc_df <- left_join(pc_df, covariates, by = "Sample") |> 
  mutate(Condition = ifelse(Condition == "CTL", "PBS", "FNF"))
pc_df$Sex <- factor(pc_df$Sex, levels = c("M", "F"))
pc_df$Condition <- factor(pc_df$Condition, levels = c("PBS", "FNF"))

expression_pca <- ggplot(data = pc_df, aes(x = PC1, 
                         y = PC2, fill = Condition, pch = Sex)) +
  geom_point(size = 2, color = "black") +
  scale_x_continuous(breaks = seq(-60, 60, 20),
                     name = paste0("PC1: ", 
                                   round(percentVar[1]*100, digits = 1), 
                                   "% variance explained")) +
  scale_y_continuous(limits = c(-60, 60), breaks = seq(-60, 60, 20),
                     name = paste0("PC2: ", 
                                   round(percentVar[2]*100, digits = 1), 
                                   "% variance explained")) +
  scale_shape_manual(values = c(21, 24)) +
  scale_fill_manual(values = c(log2fcColors[["-"]], log2fcColors[["+"]])) +
  theme_custom_scatterplot() +
  theme(legend.box = "horizontal",
        panel.background = element_rect(fill = "transparent", color = "transparent"),
        plot.background = element_rect(fill = "transparent", color = "transparent"),
        axis.text = element_text(family = "Helvetica", color = "black"),
        axis.title = element_text(family = "Helvetica", size = 10),
        legend.position = "bottom",
        legend.spacing.y = unit(0, "mm"),
        legend.box.margin = margin(-20, 0, 0, 0)) +
  guides(fill = guide_legend(label.position = "bottom", order = 1,
                             override.aes=list(shape = 21)),
         pch = guide_legend(label.position = "bottom")) +
  labs(fill = NULL, pch = NULL)

save(expression_pca, file = "plots/conditionDE_Fig1_supp/expression_pca.rda")

# Replicate correlation --------------------------------------------------------

## Replicate dataset
load("data/replicates_gse.rda")

# Put into DESeqDataSet for normalization
dds_replicate <- DESeqDataSet(gse, ~1)

# Filter out lowly expressed genes
keep <- rowSums(counts(dds_replicate) >= 10) >= ceiling(nrow(colData(gse))*0.10)
dds_replicate <- dds_replicate[keep,]

# Normalized counts
dds_replicate_norm <- assay(vst(dds_replicate))

# Split into CTL and FNF
ctl_replicate_norm <- as.data.frame(dds_replicate_norm) |> 
  dplyr::select(contains("_CTL"))
fnf_replicate_norm <- as.data.frame(dds_replicate_norm) |> 
  dplyr::select(contains("_FNF"))

## Standard dataset
# Load dds
load("data/condition_de/differential_expression_dds.rda")

# Normalized counts
dds_standard_norm <- assay(vst(dds))

# Split into CTL and FNF
ctl_standard_norm <- as.data.frame(dds_standard_norm) |> 
  dplyr::select(contains("_CTL"))
fnf_standard_norm <- as.data.frame(dds_standard_norm) |> 
  dplyr::select(contains("_FNF"))

## Calculate Pearson's correlation within donor replicates
ctl_cor_amongDonors <- cor_amongDonors(ctl_replicate_norm, "PBS")
fnf_cor_amongDonors <- cor_amongDonors(fnf_replicate_norm, "FN-f")

## Calculate Pearson's correlations across donors 
ctl_cor_acrossDonors <- cor_acrossDonors(ctl_replicate_norm,
                                         ctl_standard_norm, 
                                         "PBS")
fnf_cor_acrossDonors <- cor_acrossDonors(fnf_replicate_norm,
                                         fnf_standard_norm, 
                                         "FN-f")

# Unpaired, two-sided Wilcox test between within donor and across donor using Fisher's Z
ctl_rep_pval <- wilcox.test(ctl_cor_amongDonors$fisherz,
                            ctl_cor_acrossDonors$fisherz,
                            alternative = "two.sided",
                            paired = FALSE)$p.value
fnf_rep_pval <- t.test(fnf_cor_amongDonors$fisherz,
                       fnf_cor_acrossDonors$fisherz,
                       alternative = "two.sided",
                       paired = FALSE)$p.value


rep_pvals <- data.frame("Condition" = c("PBS", "FN-f"),
                        "pval" = c(ctl_rep_pval, fnf_rep_pval)) |> 
  mutate(label = paste0("pval = ", round(pval, digits = 3))) |> 
  mutate(Condition = factor(Condition, levels = c("PBS", "FN-f")))

# Join all correlations into one for plotting
replicate_correlations <- bind_rows(ctl_cor_amongDonors,
                                    fnf_cor_amongDonors,
                                    ctl_cor_acrossDonors,
                                    fnf_cor_acrossDonors)
replicate_correlations$Condition <- factor(replicate_correlations$Condition,
                                           levels = c("PBS", "FN-f"))


# Get n of each group
groupCounts <- replicate_correlations |> 
  group_by(Condition, group) |> 
  summarise(num = dplyr::n())

groupCountlabels <- paste0(groupCounts$group, " donors\n", "(n = ", groupCounts$num, ")")

replicate_correlation <- ggplot(replicate_correlations, aes(x = group, y = coef, fill = Condition)) +
  geom_violin() +
  stat_boxplot(geom = "errorbar", width = 0.2, linewidth = 0.75) +
  facet_wrap(~Condition, scales = "free_x") +
  scale_fill_manual(values = c(log2fcColors[["-"]], log2fcColors[["+"]])) +
  scale_y_continuous(name = "Correlation coefficient",
                     limits = c(0.94, 1),
                     breaks = seq(0.94, 1, 0.01)) +
  ggh4x::facetted_pos_scales(list(Condition == "PBS" ~ scale_x_discrete(labels = groupCountlabels[1:2]),
                                  Condition == "FN-f" ~ scale_x_discrete(labels = groupCountlabels[3:4]))) +
  geom_text(data = rep_pvals, aes(x = 1.5, y = 1, label = label),
            family = "Helvetica",
            size = 3.5,
            inherit.aes = FALSE) +
  theme_custom_scatterplot() +
  theme(legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.title.x = element_blank(),
        axis.text = element_text(color = "black"),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"))

save(replicate_correlation, 
     file = "plots/conditionDE_Fig1_supp/replicate_correlation.rda")


# Geno ancestry PCA -------------------------------------------------------
# get final donor list, eliminating geno-contaminated donor
donors <- read_csv("data/donorSamplesheet.csv", col_select = "Donor") |> 
  filter(Donor != "AM7352")

# Read in 1000G panel data
panel <- fread("../qtl/data/1000G/1000G_phase3.panel", data.table = FALSE)

# Read in PCA data from Eigenstrat, grabbing sample names and first 2 PC's
pcaData <- fread("data/CQTL_COA_01_GDA8_COA2_01_COA3_01_GDA8_COA4_COA5_COA6_COA7_ref_merge.pca.evec",
                 data.table = FALSE, select = 1:3, 
                 col.names = c("sample", "PC1", "PC2")) |> 
  # Filter for final donor list and 1000G panel sample
  filter(grepl(paste(donors$Donor, collapse = "|"), sample) |
           grepl(paste(panel$sample, collapse = "|"), sample)) |> 
  # Match ref panel to PCA data based on "sample" column
  left_join(panel, by = "sample") |> 
  # Rename our population name
  mutate(super_pop = ifelse(is.na(super_pop), "Current\nstudy", super_pop)) |> 
  mutate(super_pop = factor(super_pop, 
                            levels = c("AFR", "AMR", "EAS", "EUR", "SAS", "Current\nstudy")))

# Plot PC1 vs PC2, coloring panel by super population
geno_ancestry_pca <- ggplot(pcaData |> 
         filter(super_pop != "Current\nstudy"), 
       aes(x = PC1, y = PC2)) +
  geom_point(aes(color = super_pop), size = 0.5) +
  scale_color_manual(values = c("AFR" = ancestryColors[1], 
                                "AMR" = ancestryColors[2],
                                "EAS" = ancestryColors[3],
                                "EUR" = ancestryColors[4],
                                "SAS" = ancestryColors[5],
                                "Current\nstudy" = "black")) +
  scale_y_continuous(limits = c(-0.04, 0.04)) +
  scale_x_continuous(limits = c(-0.04, 0.04)) +
  geom_point(data = pcaData |> 
               filter(super_pop ==  "Current\nstudy"), 
             aes(PC1, PC2, fill = super_pop), size = 0.5) +
  theme_custom_scatterplot() +
  theme(legend.title = element_blank(),
        axis.line = element_line(color = "black"),
        axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.text = element_text(family = "Helvetica", size = 6,
                                   margin = margin(l = -5)))
save(geno_ancestry_pca, file = "plots/conditionDE_Fig1_supp/geno_ancestry_pca.rda")

# Assemble with plotgardener ----------------------------------------------


pdf("plots/conditionDE_Fig1_supp/SupFig1.pdf", width = 8, height = 9.75)
pageCreate(width = 8, height = 9.75, showGuides = FALSE)

## A - library read depth bar plot
plotText("A", x = 0.1, y = 0.2, just = c("left", "top"), fontfamily = "Helvetica",
         fontsize = 14, fontface = "bold")
plotGG(library_read_depth_barplot, x = 0.2, y = 0.2, width = 7.75, height = 3)

## B - qc violin plots
plotText("B", x = 0.1, y = 3.25, just = c("left", "top"), fontfamily = "Helvetica",
         fontsize = 14, fontface = "bold")
plotGG(qc_violins, x = 0.2, y = 3.25, width = 3.5, height = 3)

## C - pca
plotText("C", x = 0.1, y = 6.2, just = c("left", "top"), fontfamily = "Helvetica",
         fontsize = 14, fontface = "bold")
plotGG(expression_pca, x = 0.25, y = 6.25, width = 3.5, height = 3.5)

## D - replicate correlation
plotText("D", x = 3.75, y = 3.25, just = c("left", "top"), fontfamily = "Helvetica",
         fontsize = 14, fontface = "bold")
plotGG(replicate_correlation, x = 3.9, y = 3, width = 4.1, height = 3.25)

## E - genotyping PCA
plotText("E", x = 3.75, y = 6.2, just = c("left", "top"), fontfamily = "Helvetica",
         fontsize = 14, fontface = "bold")
plotGG(geno_ancestry_pca, x = 3.9, y = 6.25, width = 4.2, height = 3.5)

dev.off()


# Assemble GO/KEGG Supp table ---------------------------------------------

upsig_go <- read_csv("tables/SupTable2A.csv")
downsig_go <- read_csv("tables/SupTable2B.csv")
upsig_kegg <- read_csv("tables/SupTable2C.csv")
downsig_kegg <- read_csv("tables/SupTable2D.csv")


ss <- gs4_create(name = "SupTable2")
write_sheet(upsig_go,
            ss, sheet = "GO Terms - Upregulated")
write_sheet(downsig_go,
            ss, sheet = "GO Terms - Downregulated")
write_sheet(upsig_kegg,
            ss, sheet = "KEGG Pathways - Upregulated")
write_sheet(downsig_kegg,
            ss, sheet = "KEGG Pathways - Downregulated")


drive_mv(file = "SupTable2", path = as_dribble("CQTL paper/Figures and Tables"))
