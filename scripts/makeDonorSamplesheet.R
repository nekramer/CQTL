#!/usr/bin/R
suppressPackageStartupMessages(library("argparse"))
suppressPackageStartupMessages(library("googlesheets4"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("purrr"))
suppressPackageStartupMessages(library("readr"))

parser <- ArgumentParser()
parser$add_argument('--subset', default = 'freeze', 
                    help = 'Character describing which subset to get donor information from. Options are "pilot", "freeze", "replicate", and "all".')
parser$add_argument('--output', default = "donorSamplesheet.csv",
                    help = 'Output file path and name.')

args <- parser$parse_args()

# Read in DNAExtractionsLibraries sheet
gs4_auth("nekramer27@gmail.com")
samplesheet <- read_sheet(ss = "https://docs.google.com/spreadsheets/d/1JwLw9D6rMqhHC9BPrZebAN40Wojo-CqbMdiIPXAzkLo/edit#gid=1699779981",
                          sheet = "RNAExtractionsLibraries",
                          col_types = "ccccdddllllTTdddTTcccccc")

donorSamplesheet <- read_sheet(ss = "https://docs.google.com/spreadsheets/d/1JwLw9D6rMqhHC9BPrZebAN40Wojo-CqbMdiIPXAzkLo/edit#gid=1699779981",
                             sheet = "Donors",
                             col_types = "ccdcddcDDtcddcccc")


if (args$subset == "freeze"){
  
  # Grab donor IDs from RNA samplesheet
  donors <- samplesheet %>% 
    filter(NYGCQC == "PASS" & !is.na(Sequencing_Directory)) %>%
    mutate(RNAshippedDate = as.Date(RNAshippedDate)) %>%
    filter(RNAshippedDate <= "2022-03-10") %>%
    pull(Donor) %>% unique()
  
  new_donorSamplesheet <- donorSamplesheet %>% 
    filter(Donor %in% donors)
  
} else if (args$subset == "replicate"){
  
  multiReps <- samplesheet %>% 
    filter(Seq_Rep == 1 & NYGCQC == "PASS" & !is.na(Sequencing_Directory)) %>% 
    group_by(Donor, Condition) %>% 
    summarise(n = n()) %>% filter(n >= 2) %>%
    pull(Donor) %>% unique()
  
  new_donorSamplesheet <- donorSamplesheet %>% 
    filter(Donor %in% multiReps)
} else if (args$subset == "pilot"){
  
  pilot <- samplesheet %>%
    filter(NYGCQC == "PASS" & !is.na(Sequencing_Directory)) %>%
    mutate(RNAshippedDate = as.Date(RNAshippedDate)) %>%
    filter(RNAshippedDate <= "2020-12-02") %>%
    pull(Donor) %>% unique()
  
  new_donorSamplesheet <- donorSamplesheet %>% 
    filter(Donor %in% pilot)
  
} else if (args$subset == "all"){
  new_donorSamplesheet <- donorSamplesheet %>% 
    filter(TARCQC == "PASS")
} else {
  
  stop("Invalid subset option.")
}

# Filter out unhelpful columns (Notes, all NA)
final_donorSamplesheet <- new_donorSamplesheet %>%
  dplyr::select(-Notes, -`Notes for eliza`) %>%
  discard(~all(is.na(.x)))
 
write_csv(final_donorSamplesheet, file = args$output)