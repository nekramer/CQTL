library(tidyverse)
library(tximeta)

coldata_fnf <- read_csv("../processing/standard/samplesheet.csv") %>%
  # Get distinct samples
  distinct(Sample, .keep_all = TRUE) %>%
  ## Rename columns
  dplyr::rename(names = Sample) %>%
  # Add donor info
  left_join(read_csv("data/donorSamplesheet.csv", col_select = c("Donor", "Sex", "Age", "Race")) %>%
              mutate(Race = replace_na(Race, "Unknown")), by = "Donor")

coldata_fnf$files <- file.path("..", "processing", "standard", "output", "quant",
                               coldata_fnf$names, "quant.sf")

coldata_ctl <- coldata_fnf %>% 
  filter(Condition == "CTL")


coldata_oa <- read_csv("../processing/OA/samplesheet.csv") %>%
  # Get distinct samples
  distinct(Sample, .keep_all = TRUE) %>%
  ## Rename columns
  dplyr::rename(names = Sample) %>%
  ## Add condition column
  mutate(Condition = "OA") %>%
  ## Add donor info
  left_join(read_csv("data/donorSamplesheet_oa.csv", col_select = c("Donor", "Sex", "Age", "Race"),
                     col_types = "ccdc") %>%
              mutate(Race = replace_na(Race, "Unknown")), by = "Donor")

coldata_oa$files <- file.path("..", "processing", "OA", "output", "quant",
                               coldata_oa$names, "quant.sf")

coldata_oa_ctl <- bind_rows(coldata_ctl, coldata_oa)
se_oa <- tximeta(coldata_oa_ctl)
gse_oa <- summarizeToGene(se_oa)
save(gse_oa, file = paste0("data/", Sys.Date(), "_gse_oa.rda"))




coldata_all <- bind_rows(coldata_fnf, coldata_oa)


se_all <- tximeta(coldata_all)
gse_oa_fnf <- summarizeToGene(se_all)

save(gse_oa_fnf, file = paste0("data/", Sys.Date(), "_gse_oa_fnf.rda"))
