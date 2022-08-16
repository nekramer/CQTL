#!/usr/bin/env Rscript

library(data.table)
library(multidplyr)
library(dplyr, warn.conflicts = FALSE)
library(parallel)

args = commandArgs(trailingOnly = TRUE)

# Function to read in and assign data to dataframe
readFiles <- function(file){
    name <- gsub("_alleleCounts.csv", "", basename(file))
    data <- fread(file, data.table = FALSE)
    assign(name, data, envir = globalenv())
}
# Function to grab a dataframe's list of SNPs (in 3rd column)
grabSNPs <- function(dfName){
    df <- get(dfName)
    snps <- df[,3]
    #return(snps)
    return(data.frame("group" = rep(dfName, length(snps)), "variant" = snps))
}

# Function to subset each dataframe for a list of SNPs and write it to a new file
subsetData <- function(dfName, snpList){
    df <- get(dfName)
    df_subset <- df[which(df[,3] %in% snpList),]
    filename <- paste0("output/", dfName, "/alleleCounts/", dfName, "_commonCounts.csv")
    write.table(df_subset, file = filename, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
}


## Read in data and get list of names of dataframes
invisible(lapply(args, readFiles))
dfNames <- unlist(lapply(args, basename))
dfNames <- unlist(lapply(dfNames, gsub, pattern = "_alleleCounts.csv", replacement = ""))

## Find common list of SNPs between all datasets
allSnps <- lapply(dfNames, grabSNPs)


cluster <- new_cluster(parallel::detectCores())
cluster_library(cluster, "dplyr")

allSnps_df <- bind_rows(allSnps) %>%
    group_by(variant) %>%
    partition(cluster) %>%
    summarize(num_groups = n()) %>%
    collect()





# allSnps_df <- do.call(rbind, allSnps) %>%
#     group_by(variant) %>%
#     partition(cluster) %>%
#     summarize(num_groups = n()) %>%
#     ungroup()

write.table(allSnps_df, "alleleCount_stats.txt", quote = FALSE, row.names = FALSE, col.names = TRUE)










#commonSnps <- Reduce(intersect, allSnps)
## Write common list to a txt file
#write.table(c(commonSnps), "output/vcf/commonSnps.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

## Subset each dataframe for common list of snps and write new files
#invisible(lapply(dfNames, subsetData, snpList = commonSnps))