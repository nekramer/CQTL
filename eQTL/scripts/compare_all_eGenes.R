library(tidyverse)
library(tools)
conditions <- c("CTL", "FNF")

maxCond <- list()
for (condition in conditions){
  covariateOptions <- c(paste0("output/summary/", condition, "_none.txt"),
                        paste0("output/summary/", condition, "_genoPC_batch.txt"),
                        paste0("output/summary/", condition, "_genoPC.txt"),
                        paste0("output/summary/", condition, "_batch.txt"))
  
  maxEach <- list()
  for (option in covariateOptions){
    
    optionName <- gsub(paste0(condition, "_"), "", 
                       gsub(".txt", "", 
                            file_path_sans_ext(basename(option))))
    
    max_option <- read_csv(option) %>%
      filter(qval == max(qval) & FDR == max(FDR)) %>%
      mutate(Option = paste0(optionName, "_", Covariate_Options)) %>%
      dplyr::select(-Covariate_Options) %>%
      relocate(Option)
    
    maxEach[[optionName]] <- max_option
    
  }
  max_condition_option <- bind_rows(maxEach) %>%
    filter(qval == max(qval) & FDR == max(FDR))
  
  
  maxCond[[condition]] <- max_condition_option
  
}
