library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
condition <- args[1]

ld_data <- list()

for (fileNum in 2:length(args)){
  
  file <- args[fileNum]
  data <- read_csv(file, col_select = c("SNP_A", "SNP_B", "R2", "rsid")) %>%
    dplyr::rename("SNP_B_rsid" = "rsid") %>%
    dplyr::rename("variantID" = "SNP_A")
    
  ld_data[[fileNum]] <- data
}


all_ld_data <- bind_rows(ld_data)
write_csv(all_ld_data, file = paste0("data/", condition, "_LDbuddies.csv"))