library(tidyverse)


lowgrade_hg19 <- read_csv("data/Steinberg_2021/processed/eQTL_LowGradeCartilage_perm_sig_lead.csv",
                           col_types = "ccddddddddcddddddddddddddccc")

lowgrade_hg38 <- read_delim("data/Steinberg_2021/raw/eQTL_LowGradeCartilage_hg38.bed",
                             col_types = c("cddc"), 
                             col_names = c("chrom", "hg38pos", "hg38pos2", "genotype_id")) |> 
  dplyr::select(genotype_id, hg38pos) |> 
  right_join(lowgrade_hg19, by = "genotype_id", relationship = "many-to-many")


write_csv(lowgrade_hg38, file = "data/Steinberg_2021/processed/eQTL_LowGradeCartilage_perm_sig_lead_hg38.csv")



highgrade_hg19 <- read_csv("data/Steinberg_2021/processed/eQTL_HighGradeCartilage_perm_sig_lead.csv",
                           col_types = "ccddddddddcddddddddddddddccc")

highgrade_hg38 <- read_delim("data/Steinberg_2021/raw/eQTL_HighGradeCartilage_hg38.bed",
                             col_types = c("cddc"), 
                             col_names = c("chrom", "hg38pos", "hg38pos2", "genotype_id")) |> 
  dplyr::select(genotype_id, hg38pos) |> 
  right_join(highgrade_hg19, by = "genotype_id", relationship = "many-to-many")


write_csv(highgrade_hg38, file = "data/Steinberg_2021/processed/eQTL_HighGradeCartilage_perm_sig_lead_hg38.csv")
