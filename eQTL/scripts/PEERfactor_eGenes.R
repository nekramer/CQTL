library(tidyverse)
source("../plotting_utils.R")
library(qvalue)

resultPath <- "/work/users/n/e/nekramer/CQTL_eQTL/output/qtl/"

# Get list of possible PEER numbers in results ----------------------------
CTL_resultFiles <- list.files(resultPath, 
                              pattern = paste0("^CTL_PEER_k.*_perm1Mb_FDR_rsids\\.csv$"))

FNF_resultFiles <- list.files(resultPath, 
                              pattern = paste0("^FNF_PEER_k.*_perm1Mb_FDR_rsids\\.csv$"))

# Get number of significant eGenes in PBS and FN-f ----------------------------
peer_sig_eGenes <- list()

for (file in CTL_resultFiles){
  
  Nk <- as.numeric(gsub("k", "", unlist(str_split(file, "_"))[3]))
  RNAbatch <- "RNAKitBatch" %in% unlist(str_split(file, "_"))
  DNAbatch <- "DNAKitBatch" %in% unlist(str_split(file, "_"))
  sequencingBatch <- "RNASequencingBatch" %in% unlist(str_split(file, "_"))
  genotypingBatch <- "genoBatch" %in% unlist(str_split(file, "_"))
  
  batchGroup <- paste(c("RNAKitBatch", 
                  "DNAKitBatch", 
                  "RNAsequencingBatch", 
                  "genoBatch")[c(RNAbatch,
                                 DNAbatch,
                                 sequencingBatch,
                                 genotypingBatch)], collapse = "_")
  
  qval_num_eGenes <- read_csv(paste0(resultPath, file)) |> 
    mutate(qval = qvalue(adj_beta_pval, pi0.method = "bootstrap")$qvalue) |> 
    filter(qval < 0.05) |>  
    nrow()
  
  FDR_num_eGenes <- read_csv(paste0(resultPath, file))|> 
    filter(FDR < 0.05) |> 
    nrow()
  
  peer_sig_eGenes[[paste0(file, "FDR")]] <- data.frame("Condition" = "PBS",
                                                        "PEER" = Nk,
                                                        "eGenes" = FDR_num_eGenes,
                                                        "batchGroup" = batchGroup,
                                                        "correction" = "FDR")
}

for (file in FNF_resultFiles){
  
  Nk <- as.numeric(gsub("k", "", unlist(str_split(file, "_"))[3]))
  RNAbatch <- "RNAKitBatch" %in% unlist(str_split(file, "_"))
  DNAbatch <- "DNAKitBatch" %in% unlist(str_split(file, "_"))
  sequencingBatch <- "RNASequencingBatch" %in% unlist(str_split(file, "_"))
  genotypingBatch <- "genoBatch" %in% unlist(str_split(file, "_"))
  
  batchGroup <- paste(c("RNAKitBatch", 
                        "DNAKitBatch", 
                        "RNAsequencingBatch", 
                        "genoBatch")[c(RNAbatch,
                                       DNAbatch,
                                       sequencingBatch,
                                       genotypingBatch)], collapse = "_")
  
  qval_num_eGenes <- read_csv(paste0(resultPath, file)) |>
    mutate(qval = qvalue(adj_beta_pval, pi0.method = "bootstrap")$qvalue) |> 
    filter(qval < 0.05) |> 
    nrow()
  FDR_num_eGenes <- read_csv(paste0(resultPath, file)) |> 
    filter(FDR < 0.05) |> 
    nrow()
  
  peer_sig_eGenes[[paste0(file, "FDR")]] <- data.frame("Condition" = "FN-f",
                                                        "PEER" = Nk,
                                                        "eGenes" = FDR_num_eGenes,
                                                        "batchGroup" = batchGroup,
                                                       "correction" = "FDR")
}

peer_sig_eGenes <- peer_sig_eGenes |> 
  bind_rows() |> 
  mutate(batchGroup = ifelse(batchGroup == "", "No batches", batchGroup)) |> 
  mutate(batchGroup = ifelse(batchGroup == "RNAKitBatch_DNAKitBatch_RNAsequencingBatch_genoBatch", 
                             "RNAKitBatch_DNAKitBatch\nRNAsequencingBatch_genoBatch",batchGroup)) |> 
  mutate(batchGroup = factor(batchGroup, 
                             levels = c("DNAKitBatch",
                                        "RNAKitBatch",
                                        "RNAKitBatch_DNAKitBatch",
                                        "RNAKitBatch_DNAKitBatch\nRNAsequencingBatch_genoBatch",
                                        "No batches"))) |> 
  mutate(Condition = factor(Condition, levels = c("PBS", "FN-f")))

# Plot -------------------------------------------------------------------

batchLabels_all <- peer_sig_eGenes |> 
  filter(PEER == 50) |> 
  filter(batchGroup != "No batches") |> 
  mutate(y = case_when(batchGroup == "RNAKitBatch" ~ eGenes + 15,
                       batchGroup == "RNAKitBatch_DNAKitBatch" ~ eGenes - 20,
                       batchGroup == "DNAKitBatch" ~ eGenes -10,
                       batchGroup == "RNAKitBatch_DNAKitBatch\nRNAsequencingBatch_genoBatch" ~ eGenes))

batchLabels_none <- peer_sig_eGenes |> 
  filter(PEER == 50) |> 
  filter(batchGroup == "No batches") |> 
  mutate(y = eGenes + 20)


plot_peer_sig_eGenes <- peer_sig_eGenes |> 
  filter(correction == "FDR")

ggplot(data = plot_peer_sig_eGenes, 
       aes(x = PEER, y = eGenes, color = batchGroup)) +
  geom_vline(data = filter(plot_peer_sig_eGenes, Condition == "PBS"),
             aes(xintercept = 15), color = "grey", lty = 2) +
  geom_text(data = filter(plot_peer_sig_eGenes, Condition == "PBS"),
            aes(x = 15, y = 950, label = "15"), color = "grey", hjust = 0, 
            family = "Helvetica") +
  geom_vline(data = filter(plot_peer_sig_eGenes, Condition == "FN-f"),
             aes(xintercept = 16), color = "grey", lty = 2) + 
  geom_text(data = filter(plot_peer_sig_eGenes, Condition == "FN-f"),
            aes(x = 16, y = 950, label = "16"), color = "grey", hjust = 0, 
            family = "Helvetica") +
  geom_line(lwd = 0.75) +
  coord_cartesian(clip = "off") +
  geom_text(data = batchLabels_all, aes(x = Inf, y = y, label = batchGroup),
                  family = "Helvetica", size = 3, hjust = 0, vjust = 1) +
  geom_text(data = batchLabels_none, aes(x = Inf, y = y, label = batchGroup),
            family = "Helvetica", size = 3, fontface = "bold", hjust = 0, vjust = 1) +
  facet_wrap(~Condition, scales = "free_x") +
  scale_y_continuous(limits = c(0, 1000)) +
  scale_color_manual(values = c("grey75", "grey75", "grey75", "grey75", "#136079")) +
  xlab("Number of PEER factors") +
  ylab("Number of eGenes") +
  theme_custom_scatterplot() +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 14),
        axis.line = element_line(color = "black"),
        axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        plot.margin = unit(c(1, 10, 1, 1), "lines"),
        panel.spacing = unit(10, "lines"),
        legend.position = "none") 
  
ggsave(filename = "plots/PEERfactors/PEER_covariate_eGene_comparisons.pdf", 
       units = "in", width = 10, height = 8)
