library(tidyverse)
library(org.Hs.eg.db)
source("/pine/scr/n/e/nekramer/CQTL_gitupdates/CQTL/eQTL/scripts/utils.R")

# Read in interaction results
CTL_interactions <- readRDS('/pine/scr/n/e/nekramer/iteratePEER_RNAKitBatch/output/reQTL/CTL_interaction.rds')
FNF_interactions <- readRDS('/pine/scr/n/e/nekramer/iteratePEER_RNAKitBatch/output/reQTL/FNF_interaction.rds')

# Filter for signficant interactions with pval <= 0.05
CTL_sig <- CTL_interactions %>% keep(~.$pval <= 0.05)
FNF_sig <- FNF_interactions %>% keep(~.$pval <= 0.05)

# Parse gene_ids/variant IDs
CTL_sig_eGenes <- names(CTL_sig) %>%
  str_split(":", simplify = TRUE) %>% as.data.frame()
colnames(CTL_sig_eGenes) <- c("gene_id", "chr", "pos", "A1", "A2")
FNF_sig_eGenes <- names(FNF_sig) %>%
  str_split(":", simplify = TRUE) %>% as.data.frame()
colnames(FNF_sig_eGenes) <- c("gene_id", "chr", "pos", "A1", "A2")


# Convert gene_ids to gene symbols
CTL_sig_eGenes_symbols <- getGeneSymbol(CTL_sig_eGenes$gene_id, org_db = org.Hs.eg.db) %>%
  bind_rows() %>% write_csv(file = "/pine/scr/n/e/nekramer/iteratePEER_RNAKitBatch/output/reQTL/CTL_sig_eGenes.csv")
FNF_sig_eGenes_symbols <- getGeneSymbol(FNF_sig_eGenes$gene_id, org_db = org.Hs.eg.db) %>% 
  bind_rows() %>% write_csv(file = "/pine/scr/n/e/nekramer/iteratePEER_RNAKitBatch/output/reQTL/FNF_sig_eGenes.csv")