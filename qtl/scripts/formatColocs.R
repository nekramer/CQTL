library(tidyverse)
library(data.table)
library(coloc)
source("../utils.R")
# Functions ------------------------------------------------------------------

# Function to collapse results from coloc and filter for PP4 > 0.7
filter_colocs <- function(coloc_results){
  
  # Function to filter nested lists from colocalization results based on PP4
  keep_PP4 <- function(x){
    keep(x, ~.$summary["PP.H4.abf"] >= 0.7)
  }
  
  # Function to condense down necessary coloc information into dataframe
  collapse_coloc_result <- function(x, y){
    OAtype <- unlist(str_split(y, "_"))[1]
    qtl_rsid <- unlist(str_split(y, "_"))[2]
    return(data.frame(OA = OAtype,
                      qtl_rsid = qtl_rsid,
                      qtl_variantID = x$qtl_variantID, 
                      GWAS_lead = x$GWAS_lead,
                      GWAS_lead_variantID = x$GWAS_lead_variantID,
                      eGene_id = x$eGene_id,
                      eGene_name = x$eGene_name,
                      eGene_qval = x$eGene_qval,
                      signal = x$signal,
                      PP4 = x$summary[["PP.H4.abf"]]))

  }
  
  load(coloc_results)
  filtered_coloc <- all_colocs |> 
    map( ~ purrr::compact(.)) |>  
    keep(~length(.) != 0) |> 
    map_depth(1, keep_PP4) |> 
    list_flatten() |> 
    imap(collapse_coloc_result) |> 
    bind_rows()
  
  return(filtered_coloc)
}

# Function to go through each colocalization results and grab additional
# GWAS info about the variat
gwas_info <- function(coloc_result){

  gwas_lead_info <- fread(paste0("/work/users/n/e/nekramer/External/gwas/Boer_reprocessed/",
                          coloc_result[["OA"]], "/leads/EUR_", coloc_result[["OA"]], 
                          "_leads_ld_final.csv"), data.table = FALSE) |>
    filter(rsID == coloc_result[["GWAS_lead"]]) |> 
    # Determine risk allele based on EA, NEA, and beta: 
    # Here beta is the overall estimated effect size for EA
    # if beta is positive, risk = EA
    # if beta is negative, risk = NEA
    mutate(GWAS_risk_allele = ifelse(BETA > 0, EA, NEA)) |>
    # Change sign of beta to reflect risk allele when necessary
    mutate(BETA = ifelse(BETA < 0, -1*BETA, BETA)) |> 
    # Get risk allele frequency
    mutate(GWAS_risk_allele_frequency = ifelse(GWAS_risk_allele == EA, EAF,
                                               1 - EAF)) |> 
    # Grab GWAS stats based on OA type
    dplyr::select(GWAS_risk_allele, BETA, GWAS_risk_allele_frequency,
                  ends_with(paste0("_", coloc_result[["OA"]]))) |> 
    dplyr::rename(GWAS_risk_beta = BETA,
                  GWAS_OR = paste0("OR_", coloc_result[["OA"]]),
                  GWAS_pval = paste0("Pvalue_", coloc_result[["OA"]])) |> 
    distinct()
  
  return(gwas_lead_info)
  
}

# Function to find risk allele effect on expression in our qtl dataset
risk_effect <- function(coloc_result){
  if (coloc_result[["Condition"]] == "PBS"){
    file_cond <- "CTL"
    file_peer <- 20
  } else {
    file_cond <- "FNF"
    file_peer <- 22
  }
  
  # GWAS variantID alleles are in alphabetical order
  gwas_variantID_split <- unlist(str_split(coloc_result[["GWAS_lead_variantID"]], ":"))
  gwas_variantID_v1 <- paste(gwas_variantID_split, collapse = ":")
  gwas_variantID_v2 <- paste0(gwas_variantID_split[1], ":", 
                              gwas_variantID_split[2], ":",
                              gwas_variantID_split[4], ":",
                              gwas_variantID_split[3])
  gwas_variantIDs <- c(gwas_variantID_v1, gwas_variantID_v2)
  
  risk_allele <- coloc_result[["GWAS_risk_allele"]]
  
  # Get chrom 
  chrom <- gwas_variantID_split[1]
  
  # Read in nominal QTL data and find GWAS snp
  qtl_gwas_snp <- fread(paste0("data/eqtl/qtl_nom/", file_cond, "_PEER_k", 
               file_peer, "_genoPC_allSignals_nom1Mb_MAFs_", chrom , ".csv"), 
        data.table = FALSE) |> 
    filter(variantID %in% gwas_variantIDs & 
             gene_id == coloc_result[["eGene_id"]] & 
             signal == coloc_result[["signal"]]) 
  
  # QTL beta direction towards minor allele - check if risk allele is minor allele
  # and keep or flip beta sign accordingly
  if (qtl_gwas_snp[["ma"]] == coloc_result[["GWAS_risk_allele"]]){
    risk_allele_beta <- qtl_gwas_snp[["beta"]]
  } else {
    risk_allele_beta <- -1*qtl_gwas_snp[["beta"]]
  }
  
  return(risk_allele_beta)
  
}

# Function to calculate distance between coloc'ed lead eSNPs and eGenes
coloc_snp_eGene_distance <- function(coloc_result){
  
  if (coloc_result[["Condition"]] == "PBS"){
    file_cond <- "CTL"
    file_peer <- 20
  } else {
    file_cond <- "FNF"
    file_peer <- 22
  }
  
  # Grab GWAS variant position based on variantID
  gwas_variant_pos <- as.numeric(unlist(str_split(coloc_result[["GWAS_lead_variantID"]], ":"))[2])
  
  # Get qtl_rsid snp location based on qtl results
  snp_egene_distances <- read_csv(paste0("data/eqtl/", file_cond, 
                                     "_PEER_k", file_peer ,"_genoPC_cond1Mb_topSignals_rsID.csv")) |> 
    filter(rsID == coloc_result[["qtl_rsid"]] & 
             gene_id == coloc_result[["eGene_id"]] & 
             signal == coloc_result[["signal"]]) |> 
    mutate(qtl_tss_dist = abs(gene_start - variant_start),
           gwas_tss_dist = abs(gene_start - gwas_variant_pos)) |> 
    dplyr::select(qtl_tss_dist, gwas_tss_dist)
  
  return(snp_egene_distances)
}

# Function to get eQTL nominal pvalue and betas for QTL lead snp and GWAS lead snp in both conditions
coloc_eqtl_pval_beta <- function(coloc_result){
  if (coloc_result[["Condition"]] == "PBS"){
    file_cond <- "CTL"
    file_peer <- 20
  } else {
    file_cond <- "FNF"
    file_peer <- 22
  }
  
  # Get esnp-egene info for QTL lead variant
  esnp_egene <- read_csv(paste0("data/eqtl/", file_cond, 
                                "_PEER_k", file_peer ,"_genoPC_cond1Mb_topSignals_rsID.csv")) |> 
    filter(variantID == coloc_result[["qtl_variantID"]] & 
             gene_id == coloc_result[["eGene_id"]] &
             signal == coloc_result[["signal"]])
  
  # Call utility function
  qtl_pval_betas <- get_pvals_betas(esnp_egene[["variantID"]], esnp_egene)
  
  # Construct info for GWAS lead variant
  gwas_egene <- data.frame("gene_id" = coloc_result[["eGene_id"]],
                           "signal" = coloc_result[["signal"]])
  
  # Try GWAS variant IDs with both orders of alleles
  gwas_pval_betas_v1 <- get_pvals_betas(coloc_result[["GWAS_lead_variantID"]], 
                                        gwas_egene |>
                                          mutate(variantID = coloc_result[["GWAS_lead_variantID"]])) |> 
    dplyr::rename(gwas_beta_pbs = beta_pbs,
                  gwas_beta_fnf = beta_fnf,
                  gwas_nompval_pbs = nompval_pbs,
                  gwas_nompval_fnf = nompval_fnf)
  gwas_variant_split <- unlist(str_split(coloc_result[["GWAS_lead_variantID"]], ":"))
  gwas_variantID_v2 <- paste0(gwas_variant_split[1], ":",
                              gwas_variant_split[2], ":",
                              gwas_variant_split[4], ":",
                              gwas_variant_split[3])
  gwas_pval_betas_v2 <- get_pvals_betas(gwas_variantID_v2, gwas_egene |> 
                                          mutate(variantID = gwas_variantID_v2)) |> 
    dplyr::rename(gwas_beta_pbs = beta_pbs,
                  gwas_beta_fnf = beta_fnf,
                  gwas_nompval_pbs = nompval_pbs,
                  gwas_nompval_fnf = nompval_fnf)
  
  if (nrow(gwas_pval_betas_v1) > 0){
    coloc_pval_beta <- bind_cols(qtl_pval_betas, gwas_pval_betas_v1)
  } else if (nrow(gwas_pval_betas_v2) > 0){
    coloc_pval_beta <- bind_cols(qtl_pval_betas, gwas_pval_betas_v2)
  }
  return(coloc_pval_beta)

}

# Function to get PP4 for colocalization with signal in condition that
# colocalization was not identified in
condition_coloc <- function(coloc_result){
  
  if (coloc_result[["Condition"]] == "PBS"){
    file_cond <- "CTL"
    file_peer <- 20
    file_opp_cond <- "FNF"
    file_opp_peer <- 22
  } else {
    file_cond <- "FNF"
    file_peer <- 22
    file_opp_cond <- "CTL"
    file_opp_peer <- 20
  }
  
  chrom <- unlist(str_split(coloc_result[["GWAS_lead_variantID"]], ":"))[1]
  
  # Pull signal from opposite condition nominal results
  qtl_condition_data <- read_csv(paste0("data/eqtl/qtl_nom/", file_opp_cond,
                                     "_PEER_k", file_opp_peer,
                                     "_genoPC_allSignals_nom1Mb_MAFs_", chrom, ".csv")) |> 
    filter(gene_id == coloc_result[["eGene_id"]] & signal == coloc_result[["signal"]]) |> 
    # Update variantID with alleles in alphabetical order to match GWAS variantIDs
    rowwise() |> 
    mutate(variantID_alpha = paste0(variant_chr, ":", variant_start, ":", 
                                    paste(sort(c(A1, A2)), collapse = ":"))) |> 
    # Square beta standard error to get var beta
    mutate(var_beta = beta_se^2) 
  
  # GWAS signal based on LD
  gwas_variant_LD <- read_csv(paste0("/work/users/n/e/nekramer/External/gwas/Boer_reprocessed/", 
                                     coloc_result[["OA"]], "/leads/EUR_", coloc_result[["OA"]], 
                                   "_leads_ld_final.csv"),
                            col_types = "cdddccccddddddddddddddddddddddddddddddddddddddddddccdcccd") |> 
    filter(rsID == coloc_result[["GWAS_lead"]]) |> 
    filter(ldbuddy_R2 > 0.8)
  
  gwas_sum_chrom <- fread(paste0("/work/users/n/e/nekramer/External/gwas/Boer_reprocessed/",
                                 coloc_result[["OA"]],
                                 "/summary_stats/",
                                 coloc_result[["OA"]], "_", 
                                 chrom, ".csv"), data.table = FALSE) |> 
    # filter for variants based on LD
    filter(`CHR:hg38POS` %in% gwas_variant_LD$`ldbuddy_CHR:hg38POS` & !is.na(hg38pos)) |> 
    # GWAS gives EAF, so convert to MAF
    mutate(MAF = ifelse(EAF < 0.5, 
                        EAF, 1 - EAF)) |> 
    # Determine minor allele based on MAF and EAF
    mutate(MA = ifelse(EAF == MAF, EA, NEA)) |> 
    # Flip sign of beta if NEA is not MA for consistency with QTL data beta (relative to minor allele)
    mutate(BETA_adjusted = ifelse(NEA != MA, -1*BETA, BETA)) |>
    # Add variantID column with position and alphabetical alleles to match
    # QTL variantIDs
    rowwise() |> 
    mutate(variantID = paste0("chr", `CHR:hg38POS`, ":", 
                              paste(sort(c(EA, NEA)), collapse = ":")))
  
  
  # Sample sizes
  case_control_sizes <- read_csv("/work/users/n/e/nekramer/External/gwas/Boer_reprocessed/Case_Control_sampleSizes.csv") |> 
    filter(OAsubtype == coloc_result[["OA"]])
  fraction_cases <- case_control_sizes$Max_Cases/(case_control_sizes$Max_Cases + 
                                                    case_control_sizes$Max_Controls)
  
  gwasN <- case_control_sizes$Max_Cases + case_control_sizes$Max_Controls
  qtlN <- 101
  
  coloc_result <- tryCatch({coloc.abf(dataset1 = list(pvalues = qtl_condition_data$nom_pval,
                                                      N = qtlN,
                                                      MAF = qtl_condition_data$maf,
                                                      type = "quant",
                                                      beta = qtl_condition_data$beta,
                                                      varbeta = qtl_condition_data$var_beta,
                                                      snp = qtl_condition_data$variantID_alpha),
                                      dataset2 = list(beta = gwas_sum_chrom$BETA_adjusted,
                                                      s = fraction_cases,
                                                      N = gwasN,
                                                      type = "cc",
                                                      MAF = gwas_sum_chrom$MAF,
                                                      pvalues = gwas_sum_chrom$p,
                                                      snp = gwas_sum_chrom$variantID))},
                           error = function(cond) {NULL})
  
  if (!is.null(coloc_result)){
    pp4 <- coloc_result$summary[["PP.H4.abf"]]
  } else {
    pp4 <- NULL
  }
  
  return(pp4)
}


# Filtering and collapsing -----------------------------------------------------

colocs_250kb <- bind_rows(filter_colocs("data/colocalization/raw/CTL_PEER_k20_genoPC_allSignals_EUR_GWAS_250kb.rda") |> 
                            mutate(Condition = "PBS"),
                          filter_colocs("data/colocalization/raw/FNF_PEER_k22_genoPC_allSignals_EUR_GWAS_250kb.rda") |> 
                            mutate(Condition = "FN-f")) |> 
  filter(eGene_name != "RAD9A")

colocs_500kb <- bind_rows(filter_colocs("data/colocalization/raw/CTL_PEER_k20_genoPC_allSignals_EUR_GWAS_500kb.rda") |> 
                            mutate(Condition = "PBS"),
                          filter_colocs("data/colocalization/raw/FNF_PEER_k22_genoPC_allSignals_EUR_GWAS_500kb.rda") |> 
                            mutate(Condition = "FN-f"))

colocs_1Mb <- bind_rows(filter_colocs("data/colocalization/raw/CTL_PEER_k20_genoPC_allSignals_EUR_GWAS_1Mb.rda") |> 
                            mutate(Condition = "PBS"),
                          filter_colocs("data/colocalization/raw/FNF_PEER_k22_genoPC_allSignals_EUR_GWAS_1Mb.rda") |> 
                            mutate(Condition = "FN-f"))

colocs_06LD <- bind_rows(filter_colocs("data/colocalization/raw/CTL_PEER_k20_genoPC_allSignals_EUR_GWAS_06LD.rda") |> 
                           mutate(Condition = "PBS"),
                         filter_colocs("data/colocalization/raw/FNF_PEER_k22_genoPC_allSignals_EUR_GWAS_06LD.rda") |> 
                           mutate(Condition = "FN-f"))

colocs_08LD <- bind_rows(filter_colocs("data/colocalization/raw/CTL_PEER_k20_genoPC_allSignals_EUR_GWAS_08LD.rda") |> 
                           mutate(Condition = "PBS"),
                         filter_colocs("data/colocalization/raw/FNF_PEER_k22_genoPC_allSignals_EUR_GWAS_08LD.rda") |> 
                           mutate(Condition = "FN-f"))

# Get additional GWAS information
colocs_250kb <- colocs_250kb |> 
  bind_cols(bind_rows(apply(colocs_250kb, 1, gwas_info)))

# Find direction of effect of risk allele in QTL data
colocs_250kb$riskAllele_eqtl_beta <- apply(colocs_250kb, 1, risk_effect)


# Get nominal pvalues and betas of qtl and GWAS leads in QTL data
colocs_250kb <- colocs_250kb |> 
  bind_cols(bind_rows(apply(colocs_250kb, 1, coloc_eqtl_pval_beta)))

# Calculate lead QTL distance to eGene TSS
colocs_250kb <- colocs_250kb |> 
  bind_cols(bind_rows(apply(colocs_250kb, 1, coloc_snp_eGene_distance)))
colocs_250kb <- colocs_250kb |> 
  group_by(eGene_id) |> 
  # collapse duplicates
  mutate(qtl_tss_dist = paste(qtl_tss_dist, collapse = "\n"))

# Get corresponding PP4 testing lead/eGene signal in condition in which it was not
# colocalized. 
# Might not have been tested because:
# * eGene wasn't sig in that opposite condition: TGFA, PIK3R1, ABCA10, ALDH1A2, RNF144B, SMAD3, SLC44A2
#       - Note: some of these still have signals there, they're just not as strong and weren't called eGenes for that cond
# * The eGene had a different lead that wasn't in LD with the GWAS lead: FNF TRIOBP signal 0, FNF MAP2K6, FNF ABCA5, FNF ABCA9, FNF CARF
# Or it was tested and was filtered out because PP4 wasn't > 0.7: FNF TRIOBP signal 1, FNF USP8, PBS PSMG1 
# - Note: USP8 was right on the cutoff at 0.68 and PSMG1 was right on the cutoff at 0.69
colocs_250kb$PP4_other_cond <- apply(colocs_250kb, 1, condition_coloc)

# Restructure PP4 and PP4_other_cond to include condition in column name
colocs_250kb <- colocs_250kb |> 
  mutate(PP4_pbs = ifelse(Condition == "PBS", PP4, PP4_other_cond),
         PP4_fnf = ifelse(Condition == "FN-f", PP4, PP4_other_cond)) |> 
  dplyr::select(-PP4, -PP4_other_cond)

write_csv(colocs_250kb, file = "data/colocalization/processed/PBS_FNF_allSignals_EUR_GWAS_250Kb_colocs.csv")
