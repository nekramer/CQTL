library(tidyverse)

percentOverlap <- function(a, b){
  intersection <- length(intersect(a,b))
  percentOverlap <- intersection/length(a)
  return(percentOverlap)
}

conditions <- c("CTL", "FNF")

for (condition in conditions){
  
  donors73_egenes <- read_csv(paste0("output/qtl/", 
                               condition, 
                               "73donors_PEER_k10_genoPC_RNAKitBatch_perm1Mb_FDR.txt")) %>%
    filter(FDR <= 0.05) %>% nrow()
  print(donors73_egenes)
  freeze_egenes <- read_csv(paste0("/pine/scr/n/e/nekramer/iteratePEER_RNAKitBatch/output/qtl/",
                                   condition,
                                   "_PEER_k10_genoPC_RNAKitBatch_perm1Mb_FDR.txt")) %>%
    filter(FDR <= 0.05) %>% nrow()
  print(freeze_egenes)
}


percentOverlaps <- list()

for (condition in conditions){
  
  donors73_egenes_FDR <- read_csv(paste0("output/qtl/", 
                                         condition, 
                                         "73donors_PEER_k10_genoPC_RNAKitBatch_perm1Mb_FDR.txt")) %>%
    filter(FDR <= 0.05) %>%
    pull(gene_id)
  
  donors73_egenes_qval <- read_csv(paste0("output/qtl/", 
                                          condition, 
                                          "73donors_PEER_k10_genoPC_RNAKitBatch_perm1Mb_FDR.txt")) %>%
    filter(qval <= 0.05) %>%
    pull(gene_id)
  
  freeze_egenes_FDR <- read_csv(paste0("/pine/scr/n/e/nekramer/iteratePEER_RNAKitBatch/output/qtl/",
                                       condition,
                                       "_PEER_k10_genoPC_RNAKitBatch_perm1Mb_FDR.txt")) %>%
    filter(FDR <= 0.05) %>%
    pull(gene_id)
  
  freeze_egenes_qval <- read_csv(paste0("/pine/scr/n/e/nekramer/iteratePEER_RNAKitBatch/output/qtl/",
                                        condition,
                                        "_PEER_k10_genoPC_RNAKitBatch_perm1Mb_FDR.txt")) %>%
    filter(qval <= 0.05) %>%
    pull(gene_id)
  
  #percOverlap_FDR <- percentOverlap(freeze_egenes_FDR, donors73_egenes_FDR)
  #percOverlap_qval <- percentOverlap(freeze_egenes_qval, donors73_egenes_qval)
  percOverlap_FDR <- percentOverlap(donors73_egenes_FDR, freeze_egenes_FDR)
  percOverlap_qval <- percentOverlap(donors73_egenes_qval, freeze_egenes_qval)
  percentOverlaps[[paste0(condition, "_FDR")]] <- percOverlap_FDR
  percentOverlaps[[paste0(condition, "_qval")]] <- percOverlap_qval
  
}





for (option in CTL_options$Name){
    # Read in eGenes for set of covariate options
    eGenes_FDR <- read_csv(paste0("output/qtl/CTL_", option, "_perm1Mb_FDR.txt")) %>%
      filter(FDR <= 0.05) %>%
      pull(gene_id)

    eGenes_qval <- read_csv(paste0("output/qtl/CTL_", option, "_perm1Mb_FDR.txt")) %>%
      filter(qval <= 0.05) %>%
      pull(gene_id)

    # Calculate percent overlap
    percOverlap_FDR <- percentOverlap(eGenes_FDR, tissue_eGenes)
    percOverlap_qval <- percentOverlap(eGenes_qval, tissue_eGenes)
    CTL_percentOverlap[[paste0(option, "_FDR")]] <- percOverlap_FDR
    CTL_percentOverlap[[paste0(option, "_qval")]] <- percOverlap_qval

  }