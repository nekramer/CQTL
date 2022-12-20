library(tidyverse)
library(googlesheets4)
library(corrplot)
library(lubridate)
source("scripts/utils.R")

args <- commandArgs(trailingOnly = TRUE)
condition <- args[4]
Nk <- args[5]
validPEER <- read_csv(args[6], col_names = FALSE) %>% pull(X1)


if (Nk %in% validPEER){
  
  # Read in PEER factors and samplesheets -----------------------------------
  
  peer <- read_csv(args[1])
  donorInfo <- read_csv(args[2]) 
  samplesheet <- read_csv(args[3]) %>% 
    distinct(Sample, .keep_all = TRUE) %>% left_join(donorInfo, by = "Donor")
  
  
  # Separate by condition, order, and convert to numeric --------------------
  
  info <- samplesheet %>% filter(Condition == condition) %>%
    arrange(match(Sample, peer$Donor)) %>%
    mutate(across(where(is.character), as.factor)) %>%
    mutate(across(where(is.factor), as.numeric)) %>%
    mutate(across(where(is.POSIXct), as.numeric)) %>%
    dplyr::select(where(~sum(is.na(.x)) < 3))
  
  
  # Calculate Pearson's correlations and p-values ---------------------------
  cors <- correlationTests_cor(peer %>% dplyr::select(-Donor), info) %>% 
    na.omit()
  cor_p <- correlationTests_p(peer %>% dplyr::select(-Donor), info) %>% 
    na.omit()
  
  # Plot and save -----------------------------------------------------------
  if (length(cors) > 0){
    png(file = paste0("output/plots/", condition, "_PEER", Nk, "_correlation.png"))
    corrplot(cors, p.mat = cor_p, sig.level = 0.05)
    dev.off()
  } else {
      # Plot nothing
      file.create(file = paste0("output/plots/", condition, "_PEER", Nk, "_correlation.png"))
  }
  
} else {
  # Plot nothing
  file.create(file = paste0("output/plots/", condition, "_PEER", Nk, "_correlation.png"))
}






# peerP <- correlationTests_p(peer, peer) %>% 
#   as.data.frame() %>%
#   mutate(across(everything(), ~.x <= 0.05))
#   
# peerCorrelations <- correlationTests_cor(peer, peer) %>%
#   as.data.frame() %>%
# # Subset for significant correlations
#   replace(!peerP, NA)

