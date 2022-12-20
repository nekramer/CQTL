#!/usr/bin/R
library(tidyverse)

CTLcolor <- "darkorange1"
FNFcolor <- "#477CCD"

# Arguments ---------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
resultPath <- args[1]
FDRthreshold <- as.numeric(args[2])

# Get list of possible PEER numbers in results ----------------------------
CTL_resultFiles <- list.files(resultPath, 
                          pattern = paste0("^CTL_PEER_k.*_perm1Mb_FDR\\.txt$"))

FNF_resultFiles <- list.files(resultPath, 
                              pattern = paste0("^FNF_PEER_k.*_perm1Mb_FDR\\.txt$"))

# Get number of significant eGenes ----------------------------------------
peer_sig_eGenes <- list()

for (file in CTL_resultFiles){
  
  Nk <- as.numeric(gsub("k", "", unlist(str_split(file, "_"))[3]))
  
  num_eGenes <- read_csv(paste0(resultPath, file)) %>%
      filter(qval < FDRthreshold) %>% nrow()
  
  peer_sig_eGenes[[file]] <- data.frame("Condition" = "Control",
                                        "Nk" = Nk,
                                        "eGenes" = num_eGenes)
}

for (file in FNF_resultFiles){
  
  Nk <- as.numeric(gsub("k", "", unlist(str_split(file, "_"))[3]))
  
  num_eGenes <- read_csv(paste0(resultPath, file)) %>%
    filter(qval < FDRthreshold) %>% nrow()

  peer_sig_eGenes[[file]] <- data.frame("Condition" = "FN-f",
                                        "Nk" = Nk,
                                        "eGenes" = num_eGenes)
}

peer_sig_eGenes <- peer_sig_eGenes %>% bind_rows()


p <- ggplot(data = peer_sig_eGenes, mapping = aes(x = Nk, y = eGenes, color = Condition)) +
  geom_point() +
  geom_line() +
  theme_classic() +
  scale_color_manual(values = c(CTLcolor, FNFcolor)) +
  scale_x_continuous(breaks = seq(min(peer_sig_eGenes$Nk), max(peer_sig_eGenes$Nk), 5),
                   limits = c(min(peer_sig_eGenes$Nk), max(peer_sig_eGenes$Nk) + 2)) +
  xlab("Number of PEER factors") +
  ylab("Number of eGenes") + 
  theme(legend.title = element_blank(),
        axis.line = element_line(color = "grey27")) 


direct.label(p, list(dl.trans(x = x + 0.5), "last.qp", cex = 1))

ggsave(filename = "output/peer/PEER_eGene_comparisons.pdf", units = "in",
       width = 7, height = 5)
  
