library(tidyverse)
library(vcfR)
library(data.table)

# Genotyping data ---------------------------------------------------------

vcf <- vcfR2tidy(read.vcfR("data/CTLk20_FNFk22_genoPC_leadVars.vcf.gz",
                           verbose = FALSE))
variants <- vcf$fix |> 
  dplyr::select(POS, ID, REF, ALT)
geno_data <- vcf$gt |> 
  left_join(variants, by = "POS")


# High confidence PBS -----------------------------------------------------

# Read in FN-f eGenes
FNF_eGenes <- read_csv("data/eqtl/FNF_PEER_k22_genoPC_perm1Mb_sig_rsID.csv")

# Read in original response eGenes
pbs_response_eGenes <- read_csv("data/reqtl/CTL_sig01_reQTLs_PEER_k20_genoPC.csv")


# Iterate through original response and filter for at least 5 donors per geno
PBS_sig_reQTLs_geno <- c()
for (eGene in 1:nrow(pbs_response_eGenes)){
  
  # Check genotypes
  variant_id <- pbs_response_eGenes[eGene, ]$variantID
  geneID <-  pbs_response_eGenes[eGene, ]$gene_id
  
  var_geno_data <- geno_data |> 
    filter(ID == variant_id) |> 
    group_by(gt_GT_alleles) |> 
    summarise(n = dplyr::n())
  
  # At least 5 of each genotype
  if (nrow(var_geno_data) == 3 && all(var_geno_data$n >= 5)){
    PBS_sig_reQTLs_geno[[eGene]] <- pbs_response_eGenes[eGene,]
  }
  
}
PBS_sig_reQTLs_geno <- bind_rows(PBS_sig_reQTLs_geno)

# Get PBS and FN-f betas
PBS_sig_reQTLs_beta <- list()
for (chrom in 1:22){
  chrom_reQTLs <- PBS_sig_reQTLs_geno |> 
    filter(variant_chr == paste0("chr", chrom)) |> 
    dplyr::select(-beta)
  
  
  pbs_nom <- fread(paste0("data/eqtl/qtl_nom/CTL_PEER_k20_genoPC_nom1Mb_chr", 
                          chrom, ".csv"),
                   data.table = FALSE) |> 
    dplyr::select(gene_id, variantID, beta) |> 
    dplyr::rename(PBS_beta = beta)
  
  fnf_nom <- fread(paste0("data/eqtl/qtl_nom/FNF_PEER_k22_genoPC_nom1Mb_chr", 
                          chrom,".csv"),
                   data.table = FALSE) |> 
    dplyr::select(gene_id, variantID, beta) |> 
    dplyr::rename(FNF_beta = beta)
  
  chrom_reQTLs_beta <- chrom_reQTLs |> 
    left_join(pbs_nom) |> 
    left_join(fnf_nom)
  
  PBS_sig_reQTLs_beta[[chrom]] <-  chrom_reQTLs_beta
  
}
PBS_sig_reQTLs_beta <- PBS_sig_reQTLs_beta |> 
  bind_rows()


# Now filter for only in FN-f & PBS and FN-f beta have difference > 0.2
PBS_stringent_response <- PBS_sig_reQTLs_beta |> 
  filter(!gene_id %in% FNF_eGenes$gene_id) |> 
  filter(abs(FNF_beta - PBS_beta) > 0.2)

write_csv(PBS_stringent_response, 
          file = "data/reqtl/CTL_sig01_beta_donor_reQTLs_PEER_k20_genoPC.csv")

# High confidence FN-f ----------------------------------------------------

# Read in PBS eGenes
PBS_eGenes <- read_csv("data/eqtl/CTL_PEER_k20_genoPC_perm1Mb_sig_rsID.csv") 

# Read in original response eGenes
fnf_response_eGenes <- read_csv("data/reqtl/FNF_sig01_reQTLs_PEER_k22_genoPC.csv")

# Iterate through original response and filter for at least 5 donors per geno
FNF_sig_reQTLs_geno <- c()
for (eGene in 1:nrow(fnf_response_eGenes)){
  
  # Check genotypes
  variant_id <- fnf_response_eGenes[eGene, ]$variantID
  geneID <-  fnf_response_eGenes[eGene, ]$gene_id
  
  var_geno_data <- geno_data |> 
    filter(ID == variant_id) |> 
    group_by(gt_GT_alleles) |> 
    summarise(n = dplyr::n())
  
  # At least 5 of each genotype
  if (nrow(var_geno_data) == 3 && all(var_geno_data$n >= 5)){
    FNF_sig_reQTLs_geno[[eGene]] <- fnf_response_eGenes[eGene,]
  }
  
}
FNF_sig_reQTLs_geno <- bind_rows(FNF_sig_reQTLs_geno)

# Get PBS and FN-f betas
FNF_sig_reQTLs_beta <- list()
for (chrom in 1:22){
  chrom_reQTLs <- FNF_sig_reQTLs_geno |> 
    filter(variant_chr == paste0("chr", chrom)) |> 
    dplyr::select(-beta)
  
  
  pbs_nom <- fread(paste0("data/eqtl/qtl_nom/CTL_PEER_k20_genoPC_nom1Mb_chr", 
                          chrom, ".csv"),
                   data.table = FALSE) |> 
    dplyr::select(gene_id, variantID, beta) |> 
    dplyr::rename(PBS_beta = beta)
  
  fnf_nom <- fread(paste0("data/eqtl/qtl_nom/FNF_PEER_k22_genoPC_nom1Mb_chr", 
                          chrom,".csv"),
                   data.table = FALSE) |> 
    dplyr::select(gene_id, variantID, beta) |> 
    dplyr::rename(FNF_beta = beta)
  
  chrom_reQTLs_beta <- chrom_reQTLs |> 
    left_join(pbs_nom) |> 
    left_join(fnf_nom)
  
  FNF_sig_reQTLs_beta[[chrom]] <-  chrom_reQTLs_beta
  
}
FNF_sig_reQTLs_beta <- FNF_sig_reQTLs_beta |> 
  bind_rows()

# Now filter for only in FN-f & PBS and FN-f beta have difference > 0.2
FNF_stringent_response <- FNF_sig_reQTLs_beta |> 
  filter(!gene_id %in% PBS_eGenes$gene_id) |> 
  filter(abs(FNF_beta - PBS_beta) > 0.2)

write_csv(FNF_stringent_response, 
          file = "data/reqtl/FNF_sig01_beta_donor_reQTLs_PEER_k22_genoPC.csv")