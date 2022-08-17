library(googlesheets4)
library(dplyr)
library(readr)

GOOGLE_SHEET <- "https://docs.google.com/spreadsheets/d/1JwLw9D6rMqhHC9BPrZebAN40Wojo-CqbMdiIPXAzkLo/edit#gid=235376278"

rna <- read_sheet(ss = GOOGLE_SHEET, sheet = "RNA ExtractionsLibraries")
dna <- read_sheet(ss = GOOGLE_SHEET, sheet = "DNAExtractionsLibraries")
donors <- read_sheet(ss = GOOGLE_SHEET, sheet = "Donors")

# Donors we have genotyping data back for
genoDonors <- dna %>% filter(!is.na(`File Directory`)) %>% pull(Donor) %>% unlist()

# Pull these donors from RNA sample sheet
samplesheet <- rna %>% filter(Donor %in% genoDonors) %>% 
  filter(NYGCQC != "FAIL") %>%
  filter(Tech_Rep == 1)

samplesheet$Tech_Rep <- unlist(samplesheet$Tech_Rep)
samplesheet$Seq_Rep <- unlist(samplesheet$Seq_Rep)

write_csv(samplesheet, file = "../processing/samplesheet.csv")
