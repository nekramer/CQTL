library(tidyverse)
library(data.table)
library(coloc)

# Functions ------------------------------------------------------------------

# Function to collapse results from coloc and filter for any that got a PP4 result
collapse_colocs <- function(coloc_results){
  
  # Function to condense down necessary coloc information into dataframe
  collapse_coloc_result <- function(x, y){
    qtl_rsid <- y
    return(data.frame(qtl_rsid = qtl_rsid,
                      qtl_variantID = x$qtl_variantID, 
                      GWAS_lead = x$GWAS_lead,
                      GWAS_R2 = x$gwas_R2,
                      qtl_R2 = x$qtl_R2,
                      eGene_id = x$eGene_id,
                      eGene_name = x$eGene_name,
                      eGene_qval = x$eGene_qval,
                      signal = x$signal,
                      PP4 = x$summary[["PP.H4.abf"]]))
    
  }
  
  load(coloc_results)
  filtered_coloc <- as_colocs |> 
    map( ~ purrr::compact(.)) |>  
    keep(~length(.) != 0) |> 
    imap(collapse_coloc_result) |> 
    bind_rows()
  
  return(filtered_coloc)
}

# Function to go through each colocalization results and grab additional
# GWAS info about the variant
gwas_info <- function(coloc_result){
  
  gwas_lead_info <- read_csv("/proj/phanstiel_lab/External/gwas/AS/Cortes_2013/IGAS_2013_leads_hg38_ld.csv") |>
    filter(SNP == coloc_result[["GWAS_lead"]]) |>
    rowwise() |> 
    filter(variantID == ldbuddy_variantID) |> 
    ungroup() |> 
    mutate(GWAS_risk_allele = risk_allele) |>
    mutate(GWAS_risk_beta = ldbuddy_beta) |> 
    # Grab GWAS stats based on ancestry     
    dplyr::select(GWAS_risk_allele, GWAS_risk_beta, ldbuddy_p) |> 
    dplyr::rename(GWAS_pval = ldbuddy_p) |> 
    distinct()
  
  return(gwas_lead_info)
  
}

# Filtering and collapsing -----------------------------------------------------

as_colocs_250kb <- bind_rows(collapse_colocs("data/AS_CTL_colocs.rda") |> 
                               mutate(Condition = "PBS"),
                             collapse_colocs("data/AS_FNF_colocs.rda") |> 
                               mutate(Condition = "FN-f"))


# Get additional GWAS information
as_colocs_250kb <- as_colocs_250kb |> 
  bind_cols(bind_rows(apply(as_colocs_250kb, 1, gwas_info)))


write_csv(as_colocs_250kb, 
          file = "data/AS_PBS_FNF_coloc_results_250kb.csv")
