source("scripts/utils.R")
library(plyranges)
library(plotgardener)
library(data.table)

load_object <- function(file) {
  tmp <- new.env()
  load(file = file, envir = tmp)
  tmp[[ls(tmp)[1]]]
}


CTL_AI_sig <- fread("data/AIsigCTL05_final.csv", data.table = FALSE)
FNF_AI_sig <- fread("data/AIsigFNF05_final.csv", data.table = FALSE)
# CTL <- reformatRSIDdata("data/2023-01-10_AIsigCTL05_rsids.csv")
# FNF <- reformatRSIDdata("data/2023-01-10_AIsigFNF05_rsids.csv")
# 
# # DESeq results --------------------------------------------------------
# load("data/2023-01-10_AIsigCTL_.05.rda")
# load("data/2023-01-10_AIsigFNF_.05.rda")
# 
# # Reformat and join with rsid/gene info
# sigCTL_df <- reformatDESeqdata(sigCTL_.05, CTL)
# sigFNF_df <- reformatDESeqdata(sigFNF_.05, FNF)
# 
# CTL_snps_only <- sigCTL_df[which(!sigCTL_df$rsid %in% sigFNF_df$rsid),]
# FNF_snps_only <- sigFNF_df[which(!sigFNF_df$rsid %in% sigCTL_df$rsid),]
# both_snps_only <- sigCTL_df[which(sigCTL_df$rsid %in% sigFNF_df$rsid),]


# Go through and compare to Boer  -----------------------------------------

OAsubtypes <- list.files("/proj/phanstiel_lab/External/gwas/OA/Boer_2021_hg19/")
OAsubtypes <- OAsubtypes[OAsubtypes != 'README']
OAsubtypes <- OAsubtypes[OAsubtypes != 'other']

for (subtype in OAsubtypes){
  
  gwas_ld <- load_object(file.path("/proj/phanstiel_lab/External/gwas/OA/Boer_2021_hg19/", subtype, "LD",
            paste0("LD_Boer_2021_", subtype, "_rsids.rda"))) %>%
    as.list() %>%
    plyranges::bind_ranges(.id = "leadsnp")
  
  CTL_Boer_overlaps <- CTL_snps_only[which(CTL_snps_only$rsid %in% gwas_ld$rsid),] %>% 
    as.data.frame() %>%
    left_join(as.data.frame(gwas_ld), by = "rsid")
  
  FNF_Boer_overlaps <- FNF_snps_only[which(FNF_snps_only$rsid %in% gwas_ld$rsid),] %>%
    as.data.frame() %>%
    left_join(as.data.frame(gwas_ld), by = "rsid")
  
  both_Boer_overlaps <- both_snps_only[which(both_snps_only$rsid %in% gwas_ld$rsid),] %>%
    as.data.frame() %>%
    left_join(as.data.frame(gwas_ld), by = "rsid")
  
  write_csv(CTL_Boer_overlaps, file = paste0("data/", subtype, "_CTLsig_rsidOverlap.csv"))
  write_csv(FNF_Boer_overlaps, file = paste0("data/", subtype, "_FNFsig_rsidOverlap.csv"))
  write_csv(both_Boer_overlaps, file = paste0("data/", subtype, "_bothsig_rsidOverlap.csv"))
  
}


# locus zoom --------------------------------------------------------------



# OAsubtypes <- list.files("/proj/phanstiel_lab/External/gwas/OA/Boer_2021_hg19/")
# OAsubtypes <- OAsubtypes[OAsubtypes != 'README']
# OAsubtypes <- OAsubtypes[OAsubtypes != 'other']
# 
# for (subtype in OAsubtypes){
#   
#   # Read in summary stats
#   gwas_summary <- fread(file.path("/proj/phanstiel_lab/External/gwas/OA/Boer_2021_hg19/", 
#                                         subtype, 
#                                         "summary_statistics",
#                                         paste0(subtype, "_chr5.txt.gz")),
#                         data.table = FALSE)
#   
#   # Get LD buddies
#   ld_info <- load_object(file.path("/proj/phanstiel_lab/External/gwas/OA/Boer_2021_hg19/",
#                              subtype,
#                              "LD",
#                              paste0("LD_Boer_2021_", subtype, "_rsids.rda"))) %>%
#     as.list() %>% plyranges::bind_ranges(.id = "leadsnp") %>%
#     filter(seqnames == 5) %>% as.data.frame() 
#   ld_info$snp <- paste0(ld_info$seqnames, ":", ld_info$start)
#   
#   # Join
#   
#   all_data <- left_join(gwas_summary, ld_info[,c("snp", "R2", "leadsnp")], by = "snp")
#   all_data <- as.data.frame(dplyr::group_by(all_data,
#                                         LDgrp = cut(all_data$R2,
#                                                     c(0, 0.2, 0.4, 0.6, 0.8, 1))))
#   leadVar <- 170871074
#   pageCreate(height = 3.5, width = 6.5)
#   
#   
#   AI_var <- all_data %>% filter(snp == "5:171849471")
#   
#  man1 <- plotManhattan(data = all_data,
#                chrom = "chr5", 
#                #chromstart = leadVar - 10000000,
#                #chromend = leadVar + 10000000,
#                chromstart = 160871481,
#                chromend = 180789530,
#                assembly = "hg19",
#                fill = colorby("LDgrp", 
#                               palette = colorRampPalette(c("#272C73", 
#                                                            "#98CEEF", 
#                                                            "#4A9A51", 
#                                                            "#EEA740", 
#                                                            "#DD3931", 
#                                                            "grey"))),
#                leadSNP = list(snp = "5:170871074",
#                               fill = "#DD3931", fontsize = 0, pch = 18, cex = 0.75),
#                sigLine = TRUE,
#                x = 0.5, y = 0.25, width = 5, height = 1.5)
#  annoGenomeLabel(plot = man1, x = 0.5, y = 1.75)
#  
#  
#  
#  
#   
#  plotManhattan(data = all_data,
#                chrom = "chr5", 
#                chromstart = leadVar - 5000000,
#                chromend = leadVar + 5000000,
#                assembly = "hg19",
#                leadSNP = list(snp = "5:171849471",
#                               fill = "yellow", fontsize = 0, pch = 18, cex = 0.75),
#                fill = colorby("LDgrp", 
#                               palette = colorRampPalette(c("#272C73", 
#                                                            "#98CEEF", 
#                                                            "#4A9A51", 
#                                                            "#EEA740", 
#                                                            "#DD3931", 
#                                                            "grey"))),
#                x = 0.5, y = 0.25, width = 5, height = 1.5)
#  
#  
# }
