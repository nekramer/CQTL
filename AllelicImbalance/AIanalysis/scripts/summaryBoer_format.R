library(data.table)
library(stringr)
library(dplyr)
library(plyranges)
library(readr)
library(tidyr)

args = commandArgs(trailingOnly = TRUE)
# args[1]: Boer rsid data
# args[2]: Boer summary stats
# args[3]: OA subtype
# args[4]: chrom

data <- fread(args[1], data.table = FALSE, col.names = c("chrom", "pos", "rsid"))
# Remove extra characters
data[,3] <- gsub("\\[|\\]|'", "", data[,3])

dataExpanded <- data %>%
  # Split multiple rsids
  mutate(rsid = str_split(rsid, " ")) %>% 
  # Give snps with multiple rsids their own rows
  unnest(rsid)

# Read in original summary stats to join back rsids
summary_stats <- fread(args[2], data.table = FALSE)

# Join
summary_rsids <- left_join(summary_stats, dataExpanded)

# Write to file
write_csv(summary_rsids, 
          file = paste0("data/raw/", args[3], "_chr", args[4], "_summary05Final.csv"))
