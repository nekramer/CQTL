library(tidyverse)
library(qvalue)
source('/pine/scr/n/e/nekramer/CQTL_gitupdates/CQTL/eQTL/scripts/utils.R')

# Arguments ---------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
gtex_egene_path <- args[2]
gtex_signif_path <- args[3]
threshold <- args[4]
condition <- args[5]


# Read in our hits and subset for significant hits
sigData <- read_csv(args[1]) %>%
  dplyr::filter(qval < threshold) %>%
  # Add a column to make eGene:eSNP identifier
  mutate(eGene_eSNP = paste0(gene_id, ":", variantID))

# Parse GTEx --------------------------------------------------------------

gtex_tissues <- gsub(".v8.egenes.txt.gz", "", list.files(gtex_egene_path))

# eGene pi1 comparisons ----------------------------------------------------

gtex_eGene_pi1 <- list()

for (tissue in gtex_tissues){
  # Read in GTEx dataset and grab eGene overlaps
  gtex <- read_delim(paste0(gtex_egene_path, tissue, ".v8.egenes.txt.gz"),
                     delim = "\t") %>%
    # Remove decimals from ENS IDs
    mutate(across(gene_id, gsub, pattern = "\\..*", replacement = "")) %>%
    semi_join(sigData, by = "gene_id")

  # Estimate pi1 from 1 - pi0
  pi1 <- 1- pi0est(gtex$pval_beta)$pi0
  #pi1 <- 1- pi0est(gtex$pval_nominal)$pi0

  gtex_eGene_pi1[[tissue]] <- data.frame("tissue" = tissue,
                                   "pi1" = pi1)
}

gtex_eGene_pi1 %>% bind_rows() %>%
  write_csv(file = paste0("output/GTEx/",
                          condition, "_eGene_GTEx_pi1.csv"))


# eQTL pi1 comparisons -----------------------------------------------------

gtex_eQTL_pi1 <- list()
for (tissue in gtex_tissues){
  print(tissue)
  # Read in GTEx dataset
  gtex <- read_delim(paste0(gtex_signif_path, tissue, ".v8.signif_variant_gene_pairs.txt.gz"),
                     delim = "\t") %>%
    # Remove decimals from ENS IDs
    mutate(across(gene_id, gsub, pattern = "\\..*", replacement = "")) %>%
    # Clean variant_id column to match our variant IDs
    mutate(across(variant_id, gsub, pattern = "_", replacement = ":")) %>%
    mutate(across(variant_id, gsub, pattern = ":b38", replacement = "")) %>%
    # Add a column to make eGene:eSNP identifier
    mutate(eGene_eSNP = paste0(gene_id, ":", variant_id))
    
  # eGene_eSNP pairs that overlap
  overlapping <- gtex %>% semi_join(sigData, by = "eGene_eSNP")
  
  # For eGene_eSNP pairs that aren't found in gtex, sample p-values from uniform distribution
  non_overlapping <- sigData %>% anti_join(gtex, by = "eGene_eSNP") %>% pull(eGene_eSNP) %>%
    as.data.frame() %>%
    dplyr::rename("eGene_eSNP" = ".") %>%
    mutate(pval_nominal = runif(n = n()))
  
  gtex_nom_testing <- bind_rows(overlapping, non_overlapping)

  # Estimate pi1 from 1 - pi0
  # Use nominal p value here?
  pi1 <- 1- pi0est(gtex_nom_testing$pval_nominal)$pi0
  gtex_eQTL_pi1[[tissue]] <- data.frame("tissue" = tissue,
                                         "pi1" = pi1)
}

gtex_eQTL_pi1 %>% bind_rows() %>% 
  write_csv(file = paste0("output/GTEx/",
                          condition, "_eQTL_GTEx_pi1.csv"))

# GTEx downsampled eGene Percent overlaps ----------------------------------

# Get min number of significant eGenes from GTEx
minimum_sig_eGenes <- 1000000
for (tissue in gtex_tissues){

  # Get significant eGenes from GTEx tissue
  sig_tissue_eGenes <- read_delim(paste0(gtex_egene_path,
                                     tissue, ".v8.egenes.txt.gz"),
                              delim = "\t") %>% filter(qval <= 0.05) %>% nrow()
  if (sig_tissue_eGenes < minimum_sig_eGenes){
    minimum_sig_eGenes <- sig_tissue_eGenes
  }
}

gtex_percentOverlaps <- list()

for (tissue in gtex_tissues){
  # Read in significant eGenes, rank by qvalue, and downsample to minimimum number
  downsampled_tissue <- read_delim(paste0(gtex_egene_path,
                                          tissue, ".v8.egenes.txt.gz"),
                                   delim = "\t") %>% 
    filter(qval <= 0.05) %>%
    slice_min(order_by = qval, n = minimum_sig_eGenes) %>%
    # Remove decimals from ENS IDs
    mutate(across(gene_id, gsub, pattern = "\\..*", replacement = "")) %>%
    pull(gene_id)
  
  # Calculate percent overlap of our signficant eGenes with the downsampled GTEx tissue
  eGene_percOverlap <- percentOverlap(sigData$gene_id, downsampled_tissue)
  gtex_percentOverlaps[[tissue]] <- data.frame("tissue" = tissue,
                                               "percentOverlap" = eGene_percOverlap)
  
}

gtex_percentOverlaps %>% 
  bind_rows() %>% 
  write_csv(file = paste0("output/GTEx/",
                          condition,
                          "_eGene_GTEx_downsampled_percOverlap.csv"))
