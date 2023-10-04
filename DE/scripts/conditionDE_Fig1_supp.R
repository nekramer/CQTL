library(DESeq2)
library(tidyverse)
library(DescTools)
library(ggh4x)
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

pc_df <- left_join(pc_df, covariates, by = "Sample")
pc_df$Sex <- factor(pc_df$Sex, levels = c("M", "F"))

ggplot(data = pc_df, aes(x = PC1, 
                         y = PC2, color = Condition, pch = Sex)) +
  geom_point(size = 2) +
  scale_x_continuous(breaks = seq(-60, 60, 20),
                     name = paste0("PC1: ", 
                                   round(percentVar[1]*100, digits = 1), 
                                   "% variance explained")) +
  scale_y_continuous(limits = c(-60, 60), breaks = seq(-60, 60, 20),
                     name = paste0("PC2: ", 
                                   round(percentVar[2]*100, digits = 1), 
                                   "% variance explained")) +
  scale_color_manual(values = c(log2fcColors[["-"]], log2fcColors[["+"]])) +
  theme_custom_scatterplot() +
  theme(legend.box = "horizontal",
        legend.position = "bottom",
        legend.spacing.y = unit(0.5, "mm"),
        legend.box.margin = margin(-15, 0, 0, 0)) +
  guides(color = guide_legend(label.position = "bottom", order = 1),
         pch = guide_legend(label.position = "bottom")) +
  labs(color = NULL, pch = NULL)

ggsave(filename = "plots/conditionDE_Fig1_supp/expression_pca.pdf", 
       width = 6, height = 6)

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

# Unpaired, two-sided t-test between within donor and across donor using Fisher's Z
ctl_rep_pval <- t.test(ctl_cor_amongDonors$fisherz,
                       ctl_cor_acrossDonors$fisherz,
                       paired = FALSE,
                       alternative = "two.sided")$p.value
fnf_rep_pval <- t.test(fnf_cor_amongDonors$fisherz,
                       fnf_cor_acrossDonors$fisherz,
                       paired = FALSE,
                       alternative = "two.sided")$p.value


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


ggplot(replicate_correlations, aes(x = group, y = coef, fill = Condition)) +
  geom_boxplot() +
  stat_boxplot(geom = "errorbar", width = 0.2) +
  facet_wrap(~Condition, scales = "free_x") +
  scale_fill_manual(values = c(log2fcColors[["-"]], log2fcColors[["+"]])) +
  scale_y_continuous(name = "Correlation coefficient",
                     limits = c(0.94, 1),
                     breaks = seq(0.94, 1, 0.01)) +
  ggh4x::facetted_pos_scales(list(Condition == "PBS" ~ scale_x_discrete(labels = groupCountlabels[1:2]),
                                  Condition == "FN-f" ~ scale_x_discrete(labels = groupCountlabels[3:4]))) +
  geom_text(data = rep_pvals, aes(x = 1.5, y = 1, label = label),
            family = "Helvetica",
            inherit.aes = FALSE) +
  theme_custom_scatterplot() +
  theme(legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.text = element_text(color = "black"),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"))

ggsave(filename = "plots/conditionDE_Fig1_supp/replicate_correlation.pdf",
       width = 5, height = 5)


