library(tidyverse)
library(data.table)
library(coloc)
# Functions ------------------------------------------------------------------

# Function to collapse results from coloc and filter for any that got a PP4 result
collapse_colocs <- function(coloc_results){
  
  # Function to condense down necessary coloc information into dataframe
  collapse_coloc_result <- function(x, y){
    RAtype <- unlist(str_split(y, "_"))[1]
    ancestry <- unlist(str_split(y, "_"))[2]
    qtl_rsid <- unlist(str_split(y, "_"))[length(unlist(str_split(y, "_")))]
    return(data.frame(RA = RAtype,
                      ancestry = ancestry,
                      qtl_rsid = qtl_rsid,
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
  filtered_coloc <- ra_colocs |> 
    map( ~ purrr::compact(.)) |>  
    keep(~length(.) != 0) |> 
    list_flatten() |> 
    imap(collapse_coloc_result) |> 
    bind_rows()
  
  return(filtered_coloc)
}

# Function to go through each colocalization results and grab additional
# GWAS info about the variant
gwas_info <- function(coloc_result){
  
  
  if (coloc_result[["ancestry"]] == "multi"){
    ancestry <- "multi_ancestry"
  } else {
    ancestry <- coloc_result[["ancestry"]]
  }
  
  gwas_lead_info <- read_csv(paste0("/proj/phanstiel_lab/External/gwas/RA/",
                                 coloc_result[["RA"]], "/", ancestry,
                                 "/leads/", coloc_result[["RA"]], "_",
                                 ancestry, "_leads_ld_hg38_final.csv.gz")) |>
    filter(rsID == coloc_result[["GWAS_lead"]]) |>
    rowwise() |> 
    filter(variantID == ldbuddy_variantID) |> 
    ungroup() |> 
    # Determine risk allele based on EA, NEA, and beta: 
    # Here beta is the overall estimated effect size for EA
    # if beta is positive, risk = EA
    # if beta is negative, risk = NEA
    mutate(GWAS_risk_allele = ifelse(ldbuddy_beta > 0, EA, NEA)) |>
    # Change sign of beta to reflect risk allele when necessary
    mutate(GWAS_risk_beta = ifelse(ldbuddy_beta < 0, -1*ldbuddy_beta, ldbuddy_beta)) |> 
    # Get risk allele frequency
    mutate(GWAS_risk_allele_frequency = ifelse(GWAS_risk_allele == EA, as.numeric(ldbuddy_EAF),
                                               1 - as.numeric(ldbuddy_EAF))) |> 
    # Grab GWAS stats based on ancestry     
    dplyr::select(GWAS_risk_allele, GWAS_risk_beta, GWAS_risk_allele_frequency,
                  ends_with(paste0("_", coloc_result[["ancestry"]]))) |> 
    dplyr::rename(GWAS_OR = paste0("OR_", coloc_result[["ancestry"]]),
                  GWAS_pval = paste0("Pvalue_", coloc_result[["ancestry"]])) |> 
    distinct()
  
  return(gwas_lead_info)
  
}


# Filtering and collapsing -----------------------------------------------------

ra_colocs_250kb <- bind_rows(collapse_colocs("data/RA_CTL_colocs.rda") |> 
                            mutate(Condition = "PBS"),
                          collapse_colocs("data/RA_FNF_colocs.rda") |> 
                            mutate(Condition = "FN-f"))


# Get additional GWAS information
ra_colocs_250kb <- ra_colocs_250kb |> 
  bind_cols(bind_rows(apply(ra_colocs_250kb, 1, gwas_info)))


write_csv(ra_colocs_250kb, 
          file = "data/RA_PBS_FNF_coloc_results_250kb.csv")
