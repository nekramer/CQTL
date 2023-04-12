library(readr)
library(tidyr)
library(dplyr)
library(stringr)
library(purrr)
source('scripts/utils.R')

# Read in and format CTL/FNF AI results --------------------------------------

CTL <- reformatRSIDdata("data/2022-08-24_AIsigCTL_rsids.csv")
FNF <- reformatRSIDdata("data/2022-08-24_AIsigFNF_rsids.csv")

# Boer gene comparison ----------------------------------------------

Boer_S10_genes <- read_csv("data/raw/Boer_et_al_2021_S10.csv", skip = 1) %>%
  filter(!is.na(Gene)) %>% pull(Gene)

CTL_Boer_genes <- CTL[which(CTL$gene_symbol %in% Boer_S10_genes),]
FNF_Boer_genes <- FNF[which(FNF$gene_symbol %in% Boer_S10_genes),]

# Hollander comparison ----------------------------------------------------

Hollander_table2_genes <- read_csv("data/raw/Hollander_et_al_2019_Table2.csv", 
                                   skip = 1) %>%
  filter(!is.na(`Positional gene`)) %>% pull(`Positional gene`)

CTL_Hollander_genes <- CTL[which(CTL$gene_symbol %in% Hollander_table2_genes),]
FNF_Hollander_genes <- FNF[which(FNF$gene_symbol %in% Hollander_table2_genes),]


# Boer LD rsid comparison -------------------------------------------------

findCommonSNPs <- function(LDGRange, data){
  if(any(data$rsid %in% LDGRange$rsid)){
    
    # Keep snps from data
    matching_data <- data[which(data$rsid %in% LDGRange$rsid),]
    matching_LDGRange <- LDGRange[which(LDGRange$rsid %in% data$rsid),]
    
    # Merge 
    matchingAll <- merge(matching_data, 
          as.data.frame(matching_LDGRange)[,c("R2", "EffectAllele", "AlternateAllele", "p", "rsid")])
    return(matchingAll)
  } 
  
}


OA_subtypes <- list.dirs("/proj/phanstiel_lab/External/gwas/OA/Boer_2021_hg19/", 
                         recursive = FALSE, full.names = FALSE)
# Remove other
OA_subtypes <- OA_subtypes[OA_subtypes != "other"]
for (subtype in OA_subtypes){
  name <- paste0("LD_Boer_2021_", subtype, "_rsids")
  #load(paste0("/proj/phanstiel_lab/External/gwas/OA/Boer_2021_hg19/", subtype, "/LD/LD_Boer_2021_", subtype, "_rsids.rda"))
  
  #CTL_matchingSNPs <- lapply(get(name), findCommonSNPs, data = CTL)
  #CTL_matchingSNPs <- CTL_matchingSNPs %>% discard(is.null) %>% bind_rows(.id = "leadSNP")

  #FNF_matchingSNPs <- lapply(get(name), findCommonSNPs, data = FNF)
  #FNF_matchingSNPs <- FNF_matchingSNPs %>% discard(is.null) %>% bind_rows(.id = "leadSNP")
  
  #save(CTL_matchingSNPs, file = paste0("data/", subtype, "_CTL_AIoverlap.rda"))
  #save(FNF_matchingSNPs, file = paste0("data/", subtype, "_FNF_AIoverlap.rda"))
  
  load(paste0("data/", subtype, "_CTL_AIoverlap"))
  
}




























