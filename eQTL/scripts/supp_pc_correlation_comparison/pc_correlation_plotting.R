#!/usr/bin/R
library(tidyverse)
library(reshape2)
library(RColorBrewer)


# Arguments ---------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

pcs <- read_csv(args[1])
varExplained <- read_csv(args[2])
cors <- read_csv(args[3]) %>% column_to_rownames(var = "covariate")
cors_p <- read_csv(args[4]) %>% column_to_rownames(var = "covariate")
peer_cors <- read_csv(args[5]) %>% column_to_rownames(var = "covariate")
peer_cors_p <- read_csv(args[6])  %>% column_to_rownames(var = "covariate")
condition <- args[7]


# PC1 vs PC2 --------------------------------------------------------------

# Go through names of all covariates to color by
covariates <- colnames(pcs)[grep("PC", colnames(pcs), invert = TRUE)]

covariate_pc_plots <- list()
for (covariate in covariates){
  # Fix covariate name
  covariateName <- str_flatten(unlist(str_split(str_remove(covariate, "/"), " ")), "_")
  
  # Filter data for PC's and that covariate
  pcs_covariate <- pcs %>% 
    dplyr::select(contains("PC"), all_of(covariate)) %>%
    dplyr::rename("covariate" = covariate) %>%
    mutate(across(covariate, as.factor))
  
  # Plot
  p <- ggplot(data = pcs_covariate, 
              mapping = aes(x = PC1, y = PC2, color = covariate)) +
    geom_point() + 
    theme_minimal() +
    xlab(paste0("PC1 (", round(varExplained[1, "percentVarianceExplained"], digits = 2), " %)")) +
    ylab(paste0("PC2 (", round(varExplained[2, "percentVarianceExplained"], digits = 2), " %)")) +
    guides(color = guide_legend(title = covariateName)) +
    labs(subtitle = condition)
  
  covariate_pc_plots[[covariate]] <- p

}

# Save list of plots to RData for later arranging
save(covariate_pc_plots, file = paste0("output/pca/", condition, "_covarPC.RData"))

# Scree plots -------------------------------------------------------------

# Subset for first one hundred PC's
varExplained <- varExplained[1:100,]

varExplained$PC <- factor(varExplained$PC, levels = varExplained$PC)

ggplot(data = varExplained, mapping = aes(x = PC, 
                                          y = percentVarianceExplained, 
                                          group = 1)) +
  geom_line() + 
  geom_point() + 
  xlab("Principal Component") +
  ylab("Percent of Variance Explained") +
  theme_minimal() + 
  scale_x_discrete(breaks = varExplained$PC[seq(0, length(varExplained$PC), 10)])
  
ggsave(filename = paste0("output/pca/", condition, "_screeplot.pdf"))



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


ggsave(filename = paste0("output/pca/", condition, "_correlation.pdf"),
       width = 9, height = 7, units = "in")

# PEER Correlation matrix ------------------------------------------------------

covars <- c("RNAextractionKitBatch", "SequencingBatch",
            "DNAReagentBatch", "GenotypingBatch", "Sex", "Race", "Age",
            "OAGradeAvg", "CauseOfDeath", "Fragment Batch")
if (condition == "ALL"){
  covars <- c(covars, "Condition")
}


# Subset for covariates of interest
peer_cors <- peer_cors %>%
  filter(grepl("PEER", rownames(.)) | rownames(.) %in% covars) %>%
  select(contains("PEER")) %>%
  as.matrix()
  
peer_cors_p <- peer_cors_p %>%
  filter(grepl("PEER", rownames(.)) | rownames(.) %in% covars) %>%
  select(contains("PEER")) %>%
  # Add p-value stars
  mutate(across(.cols = everything(), ~ifelse(. < 0.001, "***", 
                                              ifelse(. <= 0.01, "**",
                                                     ifelse(. <= 0.1, "*", NA))))) %>%
  as.matrix()

# Melt correlation matrices and join
peer_cors_melt <- melt(peer_cors) %>%
  left_join(melt(peer_cors_p) %>% rename(star = value))

# Plot  
ggplot(peer_cors_melt, aes(Var1, Var2, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                       midpoint = 0, limit = c(-1, 1), space = "Lab", 
                       name = "Pearson\nCorrelation") +
  theme_minimal() + # minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed() + 
  geom_text(aes(Var1, Var2, label = round(value, digits = 2)), color = "black", size = 3) +
  geom_text(aes(label = star),  position = position_nudge(y = 0.15), color = "black", size = 4) + 
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank())

ggsave(filename = paste0("output/peer/", condition, "_PEERcorrelation.pdf"),
       width = 9, height = 7, units = "in")
