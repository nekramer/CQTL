#!/usr/bin/R
library(tidyverse)
library(reshape2)

# Arguments ---------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

pcs <- read_csv(args[1])

print(colnames(pcs))

varExplained <- read_csv(args[2])
batchEffectCorrection <- args[3]
cors <- read_csv(args[4]) %>% column_to_rownames(var = "covariate")
cors_p <- read_csv(args[5]) %>% column_to_rownames(var = "covariate")
condition <- args[6]

# PC1 vs PC2 --------------------------------------------------------------

# Go through batch effects to color by
batchEffect_pc_plots <- list()
for (batch in c("RNAextractionKitBatch", "SequencingBatch", "DNAReagentBatch", "GenotypingBatch", "FragmentBatch")){

  # Filter data for PC's and that batch effect
  pcs_batch <- pcs %>% 
    dplyr::select(contains("PC"), all_of(batch)) %>%
    dplyr::rename("batch" = batch) %>%
    mutate(across(batch, as.factor))
  
  # Plot
  p <- ggplot(data = pcs_batch, 
              mapping = aes(x = PC1, y = PC2, color = batch)) +
    geom_point() + 
    theme_minimal() +
    xlab(paste0("PC1 (", round(varExplained[1, "percentVarianceExplained"], digits = 2), " %)")) +
    ylab(paste0("PC2 (", round(varExplained[2, "percentVarianceExplained"], digits = 2), " %)")) +
    guides(color = guide_legend(title = batch)) +
    labs(subtitle = condition)
  
  batchEffect_pc_plots[[batch]] <- p
  
}

# Save list of plots to RData for later arranging
save(batchEffect_pc_plots, file = paste0("output/pca/", condition, "_", batchEffectCorrection, "_corrected_PCplots.RData"))

# PC Correlation matrix ------------------------------------------------------

covars <- c("RNAextractionKitBatch", "SequencingBatch",
            "DNAReagentBatch", "GenotypingBatch", "Sex", "Race", "Age",
            "OAGradeAvg", "CauseOfDeath", "Fragment Batch")
if (condition == "ALL"){
  covars <- c(covars, "Condition")
}


# Subset for first 10 PC's and for covariates of interest
cors <- cors[,1:10] %>%
  filter(rownames(.) %in% covars) %>%
  as.matrix()
cors_p <- cors_p[,1:10] %>%
  filter(rownames(.) %in% covars) %>%
  # Add p-value stars
  mutate(across(.cols = everything(), ~ifelse(. < 0.001, "***", 
                                              ifelse(. <= 0.01, "**",
                                                     ifelse(. <= 0.1, "*", NA))))) %>%
  as.matrix()

# Melt correlation matrices and join
cors_melt <- melt(cors) %>%
  left_join(melt(cors_p) %>% rename(star = value))

# Plot  
ggplot(cors_melt, aes(Var1, Var2, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                       midpoint = 0, limit = c(-1, 1), space = "Lab", 
                       name = "Pearson\nCorrelation") +
  theme_minimal() + # minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed() + 
  geom_text(aes(Var1, Var2, label = round(value, digits = 2)), color = "black", size = 4) +
  geom_text(aes(label = star),  position = position_nudge(y = 0.15), color = "black", size = 5) + 
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank())

ggsave(filename = paste0("output/pca/", condition, "_", batchEffectCorrection, "_corrected_correlation.pdf"),
       width = 9, height = 7, units = "in")
