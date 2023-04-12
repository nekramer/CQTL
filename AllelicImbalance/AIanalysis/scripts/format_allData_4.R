library(DESeq2)
library(tidyverse)
source("scripts/utils.R")

# Functions -------------------------------------------------------

getAlleleCounts <- function(data, dds){
  
  # Get variant ID that would be found in dds
  variantID <- gsub(" ", "", paste0(data[[7]], ":", data[[8]], ":", data[[9]], ":", data[[10]]))
  
  # Extract count data from DESeq dds
  countData <- plotCounts(dds = dds, gene = variantID,
             intgroup = c("condition", "donor", "allele"),
             returnData = TRUE)
  rownames(countData) <- NULL
  
  # Add rsID column for joining
  countData$rsid <- data[[16]]
  
  return(countData)
  
}

# Read in and format CTL/FNF AI rsids/genes ------------------------------------

CTL <- reformatRSIDdata("data/2023-01-10_AIresCTL_rsids.csv")
FNF <- reformatRSIDdata("data/2023-01-10_AIresFNF_rsids.csv")

# DESeq results --------------------------------------------------------

load("data/2023-01-10_AIresCTL.rda")
load("data/2023-01-10_AIresFNF.rda")
load("data/dds.rda")


# Reformat and join with rsid/gene info
resCTL_df <- reformatDESeqdata(notNA_resCTL, CTL) %>%
  dplyr::rename("CTL_log2FoldChange" = "log2FoldChange") %>%
  dplyr::rename("CTL_lfcSE" = "lfcSE") %>%
  dplyr::rename("CTL_stat" = "stat") %>%
  dplyr::rename("CTL_pvalue" = "pvalue") %>%
  dplyr::rename("CTL_padj" = "padj")
resFNF_df <- reformatDESeqdata(notNA_resFNF, FNF) %>%
  dplyr::rename("FNF_log2FoldChange" = "log2FoldChange") %>%
  dplyr::rename("FNF_lfcSE" = "lfcSE") %>%
  dplyr::rename("FNF_stat" = "stat") %>%
  dplyr::rename("FNF_pvalue" = "pvalue") %>%
  dplyr::rename("FNF_padj" = "padj")

# Join together
AIresults_df <- left_join(resCTL_df, resFNF_df)


# Get counts for each SNP -----------------------------------------------------
AIresults_alleleCounts <- bind_rows(apply(AIresults_df, 1, getAlleleCounts, dds = dds))


# Expand and join allele counts with rsid/gene info
AIresults_df_alleleCounts <- left_join(AIresults_alleleCounts, AIresults_df, by = "rsid")
AIresults_df_alleleCounts$variantID <- paste0(AIresults_df_alleleCounts$chr, 
                                           ":",
                                           AIresults_df_alleleCounts$pos,
                                           ":",
                                           AIresults_df_alleleCounts$ref,
                                           ":",
                                           AIresults_df_alleleCounts$alt) 



# Write to file -----------------------------------------------------------
write_csv(AIresults_df_alleleCounts, file = "data/AIresults_final.csv")
