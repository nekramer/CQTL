library(DESeq2)
library(tidyverse)
source("../plotting_utils.R")

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