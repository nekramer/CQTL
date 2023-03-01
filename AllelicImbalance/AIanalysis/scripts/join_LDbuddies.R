library(tidyverse)
library(data.table)


args <- commandArgs(trailingOnly = TRUE)
rsIDfile <- args[1]
ldbuddyFile <- args[2]
snp <- args[3]
condition <- args[4]

# Read in and split up rsids, remove brackets, and rename column
data <- as.data.frame(read_csv(args[1]))
data[,3] <- gsub("\\[|\\]|'", "", data[,3])
data <- data %>%
  # Split multiple rsids
  mutate(`0` = str_split(`0`, " ")) %>%
  unnest(`0`) %>%
  # Rename column
  dplyr::rename(rsid = `0`)

# Join with R2 info from original ld buddy position 
ldbuddy_positions <- fread(ldbuddyFile, data.table = FALSE)
ldbuddy_positions$CHR_B <- paste0("chr", ldbuddy_positions$CHR_B) 

data <- left_join(ldbuddy_positions, data, by = c("CHR_B", "BP_B"))


write_csv(data, file = paste0("data/", snp, "_", condition, "_ldbuddies_final.csv"))