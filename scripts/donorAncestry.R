#!/usr/bin/R
suppressPackageStartupMessages(library("argparse"))
suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("dplyr"))

parser <- ArgumentParser()
parser$add_argument('--donorSamplesheet', default = "donorSamplesheet.csv",
                    help = 'Path to donorSamplesheet to determine race breakdown.')
parser$add_argument('--output', default = "donorRace.csv",
                    help = 'Output file path and name.')

args <- parser$parse_args()

donorRace <- read_csv(args$donorSamplesheet) %>% 
  group_by(Race) %>% 
  summarise(n = n()) %>%
  arrange(desc(n))

write_csv(donorRace, file = args$output)