#!/usr/bin/R
suppressPackageStartupMessages(library("argparse"))
suppressPackageStartupMessages(library("googlesheets4"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("purrr"))
suppressPackageStartupMessages(library("readr"))

parser <- ArgumentParser()
parser$add_argument('--subset', default = 'freeze', 
                    help = 'Character describing which subset to get DNA extraction information from. Options are "pilot", "freeze", "replicate", and "all".')
parser$add_argument('--output', default = "dnaSamplesheet.csv",
                    help = 'Output file path and name.')

args <- parser$parse_args()

# Read in DNAExtractionsLibraries sheet
gs4_auth("nekramer27@gmail.com")
samplesheet <- read_sheet(ss = "https://docs.google.com/spreadsheets/d/1JwLw9D6rMqhHC9BPrZebAN40Wojo-CqbMdiIPXAzkLo/edit#gid=1699779981",
                          sheet = "RNAExtractionsLibraries",
                          col_types = "ccccdddllllTTdddTTcccccc")

dnaSamplesheet <- read_sheet(ss = "https://docs.google.com/spreadsheets/d/1JwLw9D6rMqhHC9BPrZebAN40Wojo-CqbMdiIPXAzkLo/edit#gid=1699779981",
                          sheet = "DNAExtractionsLibraries",
                          col_types = "ccdlllTddddddddTTccccc")


if (args$subset == "freeze"){
  
  # Grab donor IDs from RNA samplesheet
  donors <- samplesheet %>% 
    filter(NYGCQC == "PASS" & !is.na(Sequencing_Directory)) %>%
    mutate(RNAshippedDate = as.Date(RNAshippedDate)) %>%
    filter(RNAshippedDate <= "2022-03-10") %>%
    pull(Donor) %>% unique()
  
  new_dnaSamplesheet <- dnaSamplesheet %>% 
    filter(Donor %in% donors & GENOTYPEQC == "PASS")

} else if (args$subset == "replicate"){
  
  multiReps <- samplesheet %>% 
    filter(Seq_Rep == 1 & NYGCQC == "PASS" & !is.na(Sequencing_Directory)) %>% 
    group_by(Donor, Condition) %>% 
    summarise(n = n()) %>% filter(n >= 2) %>%
    pull(Donor) %>% unique()
  
  new_dnaSamplesheet <- dnaSamplesheet %>% 
    filter(Donor %in% multiReps & GENOTYPEQC == "PASS")
} else if (args$subset == "pilot"){
  
  pilot <- samplesheet %>%
    filter(NYGCQC == "PASS" & !is.na(Sequencing_Directory)) %>%
    mutate(RNAshippedDate = as.Date(RNAshippedDate)) %>%
    filter(RNAshippedDate <= "2020-12-02") %>%
    pull(Donor) %>% unique()
  
  new_dnaSamplesheet <- dnaSamplesheet %>% 
    filter(Donor %in% pilot & GENOTYPEQC == "PASS")
  
} else if (args$subset == "all"){
  new_dnaSamplesheet <- dnaSamplesheet %>% 
    filter(GENOTYPEQC == "PASS")
} else {
  
  stop("Invalid subset option.")
}

# Filter out unhelpful columns (all NA, Notes)
final_dnaSamplesheet <- new_dnaSamplesheet %>%
  dplyr::select(-Notes) %>%
  discard(~all(is.na(.x))) %>%
  # Rename Batch column
  dplyr::rename("GenotypingBatch" = Batch)

write_csv(final_dnaSamplesheet, file = args$output)
