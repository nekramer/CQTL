#!/usr/bin/R
suppressPackageStartupMessages(library("argparse"))
suppressPackageStartupMessages(library("googlesheets4"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("readr"))

parser <- ArgumentParser()
parser$add_argument('--subset', default = 'freeze', 
                    help = 'Character describing which subset to create samplesheet from. Options are "pilot", "freeze", and "all".')
parser$add_argument('--output', default = "samplesheet.csv",
                    help = 'Output file path and name.')

args <- parser$parse_args()

# Read in RNAExtractionsLibraries sheet
gs4_auth("nekramer27@gmail.com")
samplesheet <- read_sheet(ss = "https://docs.google.com/spreadsheets/d/1JwLw9D6rMqhHC9BPrZebAN40Wojo-CqbMdiIPXAzkLo/edit#gid=1699779981",
                          sheet = "RNAExtractionsLibraries",
                          col_types = "ccccdddllllTTdddTTccccc")

check_Sample_Reps <- function(Sample, samplesheet){
  
  donor <- Sample[['Donor']]
  condition <- Sample[['Condition']]
  # Filter for that sample's donor/condition and check other Tech Reps
  valid_samples <- samplesheet %>%
    filter(Donor == donor & Condition == condition) %>%
    filter(NYGCQC == "PASS" & !is.na(Sequencing_Directory))
  
  num_valid_samples <- valid_samples %>% distinct(Sample) %>% nrow()
  
  # If other Tech Reps are good, take first good one
  if (num_valid_samples > 1){
    
    sample_keep <- valid_samples %>% 
      slice_min(order_by = Tech_Rep)
  } else {
    sample_keep <- valid_samples
  }
  
  return(sample_keep)
  
}


# Get samples with more than 1 rep
extra_Reps <- samplesheet %>% 
  filter(Tech_Rep > 1 & NYGCQC == "PASS" & !is.na(Sequencing_Directory))

# Check if these samples have a good first Tech_Rep
keep_extra_reps <- apply(extra_Reps %>% distinct(Sample, .keep_all = TRUE), 1, check_Sample_Reps, samplesheet = samplesheet)

# Filter out extra replicates and add back in one per sample
samplesheet <- samplesheet %>%
  filter(NYGCQC == "PASS" & !is.na(Sequencing_Directory)) %>%
  # Remove extra_Reps
  anti_join(extra_Reps) %>%
  # Join back in whichever were kept from the extra reps
  bind_rows(keep_extra_reps) %>% 
  # Get distinct rows 
  distinct()

# Filter based on subset
if (args$subset == "freeze"){
  
  # Freeze at 79 donors 
  new_samplesheet <- samplesheet %>%
    mutate(RNAshippedDate = as.Date(RNAshippedDate)) %>%
    filter(RNAshippedDate <= "2022-03-10")
  
} else if (args$subset == "pilot"){
  
  # Pilot dataset of 26 donors
  new_samplesheet <- samplesheet %>%
    mutate(RNAshippedDate = as.Date(RNAshippedDate)) %>%
    filter(RNAshippedDate <= "2020-12-02")
  
} else if (args$subset == "all"){
  
  # Grab all samples, still pulling one tech rep
  new_samplesheet <- samplesheet
  
} else {
  
  stop("Invalid subset option.")
  
}

write_csv(new_samplesheet, file = args$output)
