library(tidyverse)


lowgrade <- read_csv("data/Steinberg_2021/processed/eQTL_LowGradeCartilage_perm_sig_lead.csv", 
                     col_select = c("genotype_id"), col_types = "c") |> 
  separate_wider_delim(cols = "genotype_id",
                       delim = ":",
                       names = c("chr", "pos"), cols_remove = FALSE) |> 
  mutate(chr = paste0("chr", chr),
         pos = as.numeric(pos),
         end = pos + 1) |> 
  relocate(genotype_id, .after = "end")

write_delim(lowgrade, file = "data/Steinberg_2021/raw/eQTL_LowGradeCartilage_hg19.bed",
            delim = "\t", col_names = FALSE)

highgrade <- read_csv("data/Steinberg_2021/processed/eQTL_HighGradeCartilage_perm_sig_lead.csv", 
                     col_select = c("genotype_id"), col_types = "c") |> 
  separate_wider_delim(cols = "genotype_id",
                       delim = ":",
                       names = c("chr", "pos"), cols_remove = FALSE) |> 
  mutate(chr = paste0("chr", chr),
         pos = as.numeric(pos),
         end = pos + 1) |> 
  relocate(genotype_id, .after = "end")
write_delim(highgrade, file = "data/Steinberg_2021/raw/eQTL_HighGradeCartilage_hg19.bed",
            delim = "\t", col_names = FALSE) 
