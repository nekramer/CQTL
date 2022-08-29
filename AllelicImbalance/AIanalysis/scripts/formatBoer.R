library(data.table)
library(stringr)
library(dplyr)
library(plyranges)
library(tidyr)

args = commandArgs(trailingOnly = TRUE)
# args[1]: Boer rsid data
# args[2]: Boer LD buddies R object
# args[3]: name of subtype to access the R object

mergeLDbuddies <- function(LDGranges, rsids){
  # Convert to data.frame
  LDbuddies <- as.data.frame(LDGranges)
  
  # Merge with rsid data
  rsid_LDbuddies <- merge(LDbuddies, rsids)
  
  # Convert back to GRanges
  rsidLD_GRanges <- as_granges(rsid_LDbuddies)
  
  return(rsidLD_GRanges)
}


data <- fread(args[1], data.table = FALSE, col.names = c("seqnames", "start", "rsid"))
# Remove extra characters
data[,3] <- gsub("\\[|\\]|'", "", data[,3])

dataExpanded <- data %>%
  # Remove chr from seqnames
  mutate(seqnames = str_remove(seqnames, "chr")) %>%
  # Split multiple rsids
  mutate(rsid = str_split(rsid, " ")) %>% 
  # Give snps with multiple rsids their own rows
  unnest(rsid)

# Load in original object to add back rsids
load(args[2])
objectName <- paste0("LD_Boer_2021_", args[3])

# Merge rsids with every GRanges object in the GRanges list
LD_Boer_rsids <- lapply(get(objectName), mergeLDbuddies, dataExpanded)

# Save
save(LD_Boer_rsids, 
     file = paste0("/proj/phanstiel_lab/External/gwas/OA/Boer_2021_hg19/", 
                   args[3], 
                   "/LD/LD_Boer_2021_",
                   args[3],
                   "_rsids.rda"))
