library(tidyverse)
library(janitor)
library(vcfR)
library(lme4)
source("scripts/responseQTL/testInteraction_lme4.R")

test_pair <- function(gene_variant, 
                      expression_matrix, 
                      sample_metadata, 
                      geno_matrix, 
                      qtl_formula,
                      interaction_formula){
  gene_id <- gene_variant[1]
  variant_id <- gene_variant[7]
  return(testInteractionLme4(gene_id = gene_id,
                      variant_id = variant_id,
                      expression_matrix = expression_matrix,
                      sample_metadata = sample_metadata,
                      geno_matrix = geno_matrix,
                      qtl_formula = qtl_formula, 
                      interaction_formula = interaction_formula))
}


args <- commandArgs(trailingOnly = TRUE)

# Expression --------------------------------------------------------------

# Read in combined, normalized RNA-seq counts
expression <- read_delim(args[1])

# Construct dataframe for gene meta data
gene_metadata <- expression %>%
  dplyr::select(gene_id, gene_name, `#chr`, start, end, strand)
  
expression <- expression %>%
  dplyr::select(-gene_name, -`#chr`, -start, -end, -strand)

# Construct dataframe for sample meta data
sample_metadata <- 
  data.frame("Sample" = expression %>% dplyr::select(-gene_id) %>% colnames()) %>%
  separate(col = Sample, into = c(NA, "Donor", NA, "Condition", NA, NA), 
           sep = "_",
           remove = FALSE)
  
# Convert to matrix
expression <- as.matrix(expression)
gene_ids <- expression[,"gene_id"]
expression <- expression[,-1]
expression <- apply(expression, 2, as.numeric)
rownames(expression) <- gene_ids

# Covariates --------------------------------------------------------------

# Read in covariates
covariates <- read_csv(args[2])
  
# Pull names to add to formulas
covariateNames <- colnames(covariates)[2:length(colnames(covariates))]

# Add covariates to sample_metadata
sample_metadata <- sample_metadata %>%
  left_join(covariates)

# Genotypes ---------------------------------------------------------------

# Read in filtered vcf
filtered_vcf <- read.vcfR(args[3], verbose = FALSE)
# Get into genotype matrix
genotype_matrix <- extract.gt(filtered_vcf)

# Interaction testing -----------------------------------------------------

# Define formulas
formula_qtl <- as.formula(paste("expression ~ genotype + Condition + (1|Donor) ",
                                paste(covariateNames, collapse = " + "), 
                                sep = "+ "))

formula_interaction = as.formula(paste("expression ~ genotype + Condition + Condition:genotype + (1|Donor) ",
                                       paste(covariateNames, collapse = " + "), 
                                       sep = "+ "))


# Read in eGene-variant pairs to test
eGene_variants <- read_csv(args[4])

# Test with testInteraction_lme4
interactionResults <- apply(eGene_variants, 1, test_pair,
      expression_matrix = expression,
      sample_metadata = sample_metadata,
      geno_matrix = genotype_matrix,
      qtl_formula = formula_qtl,
      interaction_formula = formula_interaction)
names(interactionResults) <- paste0(eGene_variants$gene_id, "_", 
                                    eGene_variants$variantID)


# Filter for significant hits
threshold <- as.numeric(args[5])
sig_interactionResults <- interactionResults %>% keep(~.$pval < threshold)

# Join with original QTLtools results
sig_pvals <- lapply(sig_interactionResults, `[[`, "pval")
reQTL_df <- tibble("ID" = names(sig_pvals),
                   "interaction_pval" = unlist(sig_pvals)) %>%
  separate(ID, into = c("gene_id", "variantID"), sep = "_") %>%
  left_join(eGene_variants, by = c("gene_id", "variantID"))


# Save to output
saveRDS(sig_interactionResults, file = args[6])
write_csv(reQTL_df, file = args[7])
