library(tidyverse)
library(googlesheets4)
gs4_auth("nekramer27@gmail.com")
## Genotyping/RNA-seq/QTL donors
geno_rna_samplesheet <- read_csv("qtl/data/donorSamplesheet.csv", 
                                 col_select = c("Donor", "Sex", "Age", "Race")) |> 
  mutate(Race = replace_na(Race, "Unknown")) |> 
  # Read in and join ancestries determined through genotyping pca
  left_join(read_csv("qtl/data/CQTL_COA_01_GDA8_COA2_01_COA3_01_GDA8_COA4_COA5_COA6_COA7_predictedAncestry.csv") |> 
              separate_wider_delim(cols = "Donor", delim = "_", 
                                   names = c(NA, "Donor", NA), too_many = "drop"), by = "Donor")
# Sex
n_males_rna <- geno_rna_samplesheet |> filter(Sex == "M") |> nrow()
n_females_rna <- geno_rna_samplesheet |> filter(Sex == "F") |> nrow()

# Age
age_range_rna <- range(geno_rna_samplesheet$Age)
age_mean_rna <- mean(geno_rna_samplesheet$Age)

# 1000G ancestry
perc_afr_rna <- ((geno_rna_samplesheet |> filter(Predicted_Ancestry == "AFR") |> nrow())/101) * 100
perc_amr_rna <- ((geno_rna_samplesheet |> filter(Predicted_Ancestry == "AMR") |> nrow())/101) * 100
perc_eas_rna <- ((geno_rna_samplesheet |> filter(Predicted_Ancestry == "EAS") |> nrow())/101) * 100
perc_eur_rna <- ((geno_rna_samplesheet |> filter(Predicted_Ancestry == "EUR") |> nrow())/101) * 100
perc_sas_rna <- ((geno_rna_samplesheet |> filter(Predicted_Ancestry == "SAS") |> nrow())/101) * 100

## Hi-C donors
hic_donors <- read_sheet(ss = "https://docs.google.com/spreadsheets/d/1JwLw9D6rMqhHC9BPrZebAN40Wojo-CqbMdiIPXAzkLo/edit#gid=1699779981",
                          sheet = "Hi-CLibraries",
                          col_types = "cccddddccc") |> 
  pull(Genotype) |> 
  unique()

hic_samplesheet <- read_sheet(ss = "https://docs.google.com/spreadsheets/d/1JwLw9D6rMqhHC9BPrZebAN40Wojo-CqbMdiIPXAzkLo/edit#gid=1699779981",
                              sheet = "Donors",
                              col_types = "ccdcddcDDtcddcccc") |> 
  filter(Donor %in% hic_donors) |> 
  dplyr::select(Donor, Sex, Age)

## Sex
n_males_hic <- hic_samplesheet |> filter(Sex == "M") |> nrow()
n_females_hic <- hic_samplesheet |> filter(Sex == "F") |> nrow()

## Age
age_range_hic <- range(hic_samplesheet$Age)
age_mean_hic <- mean(hic_samplesheet$Age)

## ATAC-seq donors
atac_donors <- c("AM7754", "AM7755", "AM7763")
atac_samplesheet <- read_sheet(ss = "https://docs.google.com/spreadsheets/d/1JwLw9D6rMqhHC9BPrZebAN40Wojo-CqbMdiIPXAzkLo/edit#gid=1699779981",
                              sheet = "Donors",
                              col_types = "ccdcddcDDtcddcccc") |> 
  filter(Donor %in% atac_donors) |> 
  dplyr::select(Donor, Sex, Age)

## Sex
n_males_atac <- atac_samplesheet |> filter(Sex == "M") |> nrow()
n_females_atac <- atac_samplesheet |> filter(Sex == "F") |> nrow()

## Age
age_range_atac <- range(atac_samplesheet$Age)
age_mean_atac <- mean(atac_samplesheet$Age)
