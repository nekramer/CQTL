library(tidyverse)

# Arguments ---------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

# Our dataset
eGenes <- args[1]
# Steinberg low-grade eGenes

low_grade <- read_delim(args[2], 
           delim = "\t", 
           col_types = "ccddddddddcdddddddddcc") %>% 
  distinct(Gene, .keep_all = TRUE) %>%
  separate(Gene, into = c("gene_name", "gene_id"), sep = "_")

# Steinberg high-grade eGenes
high_grade <- read_delim(args[3], 
           delim = "\t", 
           col_types = "ccddddddddcdddddddddcc") %>%
  distinct(Gene, .keep_all = TRUE) %>%
  separate(Gene, into = c("gene_name", "gene_id"), sep = "_")


#' Calculate Pi1 replicability statistic between two sets of pvalues
#' 
#' Currently expects the following columns: gene_id, qvalue and p_beta.
#' 
#' @param table1 First table maximum p-values per feature.
#' @param table2 Second table of maximum p-values per feature.
#' @param qvalue_thresh qvalue threshold for table1.
#' @return None
#' @author Kaur Alasoo
#' @export 
calculatePi1 <- function(table1, table2, qvalue_thresh = 0.1, feature_id = "gene_id"){
  #Identify significant hits from first table
  table1_hits = dplyr::filter(table1, qvalue < qvalue_thresh)
  #Extract the same genes from the second table
  table2_hits = dplyr::semi_join(table2, table1_hits, by = feature_id)
  #Estimate the proportion of replicated qtls
  pi1 = 1 - qvalue::qvalue(table2_hits$p_beta)$pi0
  return(pi1)
}

#' Calculate all pairwise Pi1 statistics for a list p-value tables.
#' 
#' Currently expects the following columns: gene_id, qvalue and p_beta.
#' 
#' @param qtl_list List of p-value tables
#' @param qvalue_thresh qvalue threshold for table1.
#' @return None
#' @author Kaur Alasoo
#' @export 
calculatePairwisePi1 <- function(qtl_list, qvalue_thresh = 0.1, tidy = FALSE, feature_id = "gene_id"){
  sample_names = names(qtl_list)
  rep_matrix = matrix(1,length(sample_names),length(sample_names))
  colnames(rep_matrix) = sample_names
  rownames(rep_matrix) = sample_names
  
  #Iterate through all pairs of p-values
  for (sn1 in 1:length(sample_names)){
    for (sn2 in 1:length(sample_names)){
      if (sn1 != sn2){
        rep_matrix[sn1, sn2] = calculatePi1(qtl_list[[sn1]], qtl_list[[sn2]], qvalue_thresh, feature_id)
      }
    }
  }
  
  #If tidy then return data frme insted of matrix
  if(tidy == TRUE){
    res = as.data.frame(rep_matrix) %>% 
      dplyr::mutate(first = rownames(rep_matrix)) %>% 
      dplyr::select(first, everything()) %>% 
      tidyr::gather("second","pi1",2:(ncol(rep_matrix)+1)) %>% 
      dplyr::arrange(first, second)
    return(res)
  }
  return(rep_matrix)
}