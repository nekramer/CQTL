library(data.table)
library(tidyverse)

OAsubtypes <- list.files("/proj/phanstiel_lab/External/gwas/OA/Boer_2021_hg19/")
OAsubtypes <- OAsubtypes[OAsubtypes != 'README']
OAsubtypes <- OAsubtypes[OAsubtypes != 'other']

for (subtype in OAsubtypes){
  
  for (chr in 1:22){
    
    chr_gwas_summary <- fread(file.path("/proj/phanstiel_lab/External/gwas/OA/Boer_2021_hg19/", 
                                    subtype, 
                                    "summary_statistics",
                                    paste0(subtype, "_chr", chr,".txt.gz")),
                          data.table = FALSE)
    
    # Grab summary stats with p-vals < 0.05 (nominal threshold)
    
    filtered_gwas_summary <- chr_gwas_summary %>%
      filter(p < 0.05)
    write_csv(filtered_gwas_summary, 
              file = paste0("data/raw/",
                            subtype, "_chr", chr, "_summary05.csv"))
    
  }
}