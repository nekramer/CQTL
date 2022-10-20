#!/usr/bin/R
library(tidyverse)
library(lubridate)
library(corrplot)
source("scripts/utils.R")

args <- commandArgs(trailingOnly = TRUE)
condition <- args[4]



# Condition-specific PCA plot ---------------------------------------------

pcs <- read_csv(args[1]) %>% filter(Donor != 'sdev')

pc_variance <- read_csv(args[1]) %>% filter(Donor == 'sdev') %>%
  dplyr::select(-Donor) %>% t() %>% as.data.frame() %>%
  rownames_to_column(var = "PC") %>%
  rename("V1" = "sdev") %>%
  mutate("PC" = as.numeric(gsub("PC", "", pc_variance$PC))) %>%
  # Calculate variance explained 
  mutate(varianceExplained = sdev^2/sum(sdev^2)) %>%
  # Percent variance explained
  mutate(percentVarianceExplained = 100*varianceExplained) 

# plot of PC1 vs PC2
ggplot(data = pcs, mapping = aes(x = PC1, y = PC2)) +
  geom_point() + theme_minimal() + 
  xlab(paste0("PC1 (", round(pc_variance$percentVarianceExplained[1], digits = 2), " %)")) +
  ylab(paste0("PC2 (", round(pc_variance$percentVarianceExplained[2], digits = 2), " %)"))

ggsave(filename = "output/plots/", condition, "_pca.pdf")


# Scree plots -------------------------------------------------------------

# Scree plot of PCs
ggplot(data = pc_variance, mapping = aes(x = PC, y = percentVarianceExplained)) +
  geom_line() + 
  geom_point() + 
  xlab("Principal Component") +
  ylab("Percent of Variance Explained") +
  scale_x_continuous(breaks = seq(0, nrow(pc_variance), 5)) + 
  theme_minimal()

ggsave(filename = "output/plots/", condition, "_screeplot.pdf")

# Pearson's correlations and p-values ------------------------------------

donorInfo <- read_csv(args[2]) 
samplesheet <- read_csv(args[3]) %>% 
  distinct(Sample, .keep_all = TRUE) %>% left_join(donorInfo, by = "Donor")

info <- samplesheet %>% filter(Condition == condition) %>%
  arrange(match(Donor, pcs$Donor)) %>%
  mutate(across(where(is.character), as.factor)) %>%
  mutate(across(where(is.factor), as.numeric)) %>%
  mutate(across(where(is.POSIXct), as.numeric)) %>%
  dplyr::select(where(~sum(is.na(.x)) < 3))

cors <- correlationTests_cor(pcs %>% dplyr::select(-Donor), info) %>% 
  na.omit()
cor_p <- correlationTests_p(pcs %>% dplyr::select(-Donor), info) %>% 
  na.omit()

# Plot and save -----------------------------------------------------------
png(file = paste0("output/plots/", condition, " _PC_correlation.png"))
corrplot(cors, p.mat = cor_p, sig.level = 0.05)
dev.off()
