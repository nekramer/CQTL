#!/usr/bin/env Rscript

library(data.table)
library(dplyr)
library(ggplot2)
library(e1071)

args = commandArgs(trailingOnly = TRUE)

## Read in pca data
pcaData <- fread(args[1], data.table = FALSE)
## Grab sample names and first 2 PC's
pcaData <- pcaData[,c(1:3)]
colnames(pcaData) <- c("sample", "PC1","PC2")

## Read in panel data
panel <- fread(args[2], data.table = FALSE)


## Match ref panel to pca data based on "sample" column
pca_panel <- left_join(pcaData, panel, by = c("sample"))

## SVN to infer super_pop
pca_panel_train <- pca_panel %>% filter(!is.na(super_pop))
pca_panel_train$super_pop <- as.factor(pca_panel_train$super_pop)
pca_panel_test <- pca_panel %>% filter(is.na(super_pop))
svm_ancestry <- svm(super_pop~PC1+PC2, data = pca_panel_train,
                    type = "C-classification", kernel = "radial")

prediction <- predict(svm_ancestry, pca_panel_test[,c("PC1", "PC2")])
pca_panel_test$super_pop <- prediction

pca_panel_test %>% dplyr::select(sample, super_pop) %>%
  dplyr::rename("Donor" = sample) %>%
  dplyr::rename("Predicted_Ancestry" = super_pop) %>%
  write_csv(file = args[5])


## Rename our population name to argument name
pca_panel[which(is.na(pca_panel$super_pop)), "super_pop"] <- args[3]

## Separate panel data from sample data to control order of plotting points
panel_df <- pca_panel[which(pca_panel$super_pop != args[3]),]
pca_df <- pca_panel[which(pca_panel$super_pop == args[3]),]

## Change super_pop columns to factors to color by values
panel_df$super_pop <- as.factor(panel_df$super_pop)
pca_df$super_pop <- as.factor(pca_df$super_pop)

## Plot PC1 vs PC2, coloring by super population
pcaplot <- ggplot(panel_df, aes(PC1, PC2)) + 
    geom_point(aes(color = factor(super_pop))) + 
    geom_point(data = pca_df, aes(PC1, PC2, fill = factor(super_pop))) +
    theme_light() + labs(color = "Population", fill = NULL)

ggsave(filename = args[4], plot = pcaplot)


