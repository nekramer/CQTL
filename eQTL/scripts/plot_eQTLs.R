#!/usr/bin/R
library(tidyverse)
library(vcfR)
# Arguments ---------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

# Response eQTL results

# Go through each response gene to compare in CTL and FNF

# Normalized expression data
CTL_normQuant <- read_delim(args[1])
FNF_normQuant <-

# Genotype data
# Read in vcf
vcf <- vcfR2tidy(read.vcfR( vcf_file, verbose = FALSE))
# Extract genotype matrix and join with info
geno_data <- left_join(vcf$gt, vcf$fix, by = c("ChromKey", "POS"))

# eQTL results (for p-value and beta value?)
  # Go through each result
  # Filter geno_matrix  for variant
  var_subset <- geno_data %>% filter(ID == "chr1:964905:C:T")
  ref <- unique(var_subset$REF)
  alt <- unique(var_subset$ALT)
  # Reformat gene_expression for that one gene for joining
  gene_expression <-  expression %>% t() %>% as.data.frame() %>% rownames_to_column(var = "Indiv")
  
  # Join with genotype data
  all_data <- left_join(var_subset, gene_expression) %>%
    # Put genotypes in factor order of ref to alt
    mutate(across(gt_GT_alleles, factor, levels = c(paste0(ref, "/", ref),
                                                    paste0(ref, "/", alt),
                                                    paste0(alt, "/", alt))))
    group_by(gt_GT_alleles) %>%
    mutate(numGeno = n())
  
  #gene_id
  #rsid
  #pval 
  #beta
  
  # Plot
  ggplot(all_data, aes(x = gt_GT_alleles, y = V1)) + 
    geom_point() +
    geom_boxplot() +
    geom_jitter(position = position_jitter(width = .1))  +
    theme_minimal() +
    ylab(label = "Normalized expression") +
    theme(axis.title.x = element_blank())
  