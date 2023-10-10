library(tidyverse)
library(data.table)
source("../plotting_utils.R")


# Read in PCA data from Eigenstrat
pcaData <- fread("data/CQTL_COA_01_GDA8_COA2_01_COA3_01_GDA8_COA4_COA5_COA6_COA7_ref_merge.pca.evec",
                 data.table = FALSE)

# Grab sample names and first 2 PC's
pcaData <- pcaData[,c(1:3)]
colnames(pcaData) <- c("sample", "PC1","PC2")

# Read in panel data
panel <- fread("data/1000G/1000G_phase3.panel", data.table = FALSE)

# Match ref panel to pca data based on "sample" column
pca_panel <- left_join(pcaData, panel, by = c("sample"))

# Rename our population name
pca_panel[which(is.na(pca_panel$super_pop)), "super_pop"] <- "Study population"

# Separate panel data from sample data to control order of plotting points
panel_df <- pca_panel[which(pca_panel$super_pop != "Study population"),] |> 
  mutate(super_pop = factor(super_pop))
pca_df <- pca_panel[which(pca_panel$super_pop == "Study population"),] |> 
  mutate(super_pop = factor(super_pop))

# Plot PC1 vs PC2, coloring panel by super population
ggplot(panel_df, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = super_pop)) +
  scale_color_manual(values = raceColors) +
  scale_y_continuous(limits = c(-0.04, 0.04)) +
  scale_x_continuous(limits = c(-0.04, 0.04)) +
  geom_point(data = pca_df, aes(PC1, PC2, fill = super_pop)) +
  theme_custom_scatterplot() +
  theme(legend.title = element_blank(),
        axis.line = element_line(color = "black"),
        axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"))

ggsave(filename = "plots/ancestry/study_ancestry.pdf",
       width = 4, height = 4)

