library(tidyverse)
library(googledrive)
library(googlesheets4)


load("data/hic/CQTL_10kb_sig_gained_loops.rda")

cqtl_10kb_diff_loops <- CQTL_10kb_sig_gained_loops |> 
  as_tibble() |> 
  dplyr::select(-width1, -strand1, -width2, -strand2, -baseMean, -pvalue) |> 
  dplyr::rename(chrom1 = seqnames1,
                chrom2 = seqnames2) |> 
  relocate(log2FoldChange, .after = end2) |> 
  relocate(lfcSE, .after = log2FoldChange) |> 
  relocate(padj, .after = lfcSE) |> 
  mutate(log2FoldChange = log2FoldChange * -1) |> 
  mutate(cluster = ifelse(log2FoldChange > 0, "gained", "lost")) |> 
  relocate(cluster, .after = padj) |> 
  rename_with(~gsub("AM7682", "Donor1", .), starts_with("AM7682")) |> 
  rename_with(~gsub("AM7683", "Donor2", .), starts_with("AM7683")) |> 
  rename_with(~gsub("AM7697", "Donor3", .), starts_with("AM7697")) |> 
  rename_with(~gsub("AM7698", "Donor4", .), starts_with("AM7698")) |> 
  rename_with(~gsub("_(\\d+)_(\\d+)_", "_r\\2_", .), starts_with("Donor")) |> 
  rename_with(~gsub("inter_30.hic", "raw", .), starts_with("Donor")) |> 
  arrange(chrom1)

write_csv(cqtl_10kb_diff_loops, "tables/SupTable12.csv")

ss <- gs4_create(name = "SupTable12")
write_sheet(cqtl_10kb_diff_loops,
            ss, sheet = "FNF_vs_PBS")


drive_mv(file = "SupTable12", path = as_dribble("CQTL paper/Figures and Tables"))
