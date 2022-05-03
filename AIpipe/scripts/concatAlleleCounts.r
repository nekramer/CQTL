#!/usr/bin/env Rscript

library(data.table)
library(stringr)
args = commandArgs(trailingOnly = TRUE)
# args[1]: donor order
# args[2]: treatment order
# args[3]...args[n]: common alleleCount files

# Convert donor/condition string inputs to vectors
donors <- unlist(str_split(args[1], ","))
conditions <- unlist(str_split(args[2], ","))

# Function to read in and assign data to dataframe
readFiles <- function(file){
    name <- gsub("_commonCounts.csv", "", basename(file))
    data <- fread(file, data.table = FALSE)
    assign(name, data, envir = globalenv())
}

# Function to reorder dataframes based on a variant vector
reorderData <- function(dfName, variantList){
    data <- get(dfName)
    data <- data[match(variantList, data$variantID),]
    assign(dfName, data, envir = globalenv())
}

# Function to find index of donor match in a vector of names
donorMatch <- function(donor, vector){
    return(str_which(vector, donor))
}

## Read in data and get list of names of dataframes
invisible(lapply(args[3:length(args)], readFiles))
dfNames <- unlist(lapply(args[3:length(args)], basename))
dfNames <- unlist(lapply(dfNames, gsub, pattern = "_commonCounts.csv", replacement = ""))

## Put all the dataframes in the same order of SNPs
# Get order from first dataframe
variantOrder <- get(dfNames[1])[,3] 
# Put other dataframes in this same order
invisible(lapply(dfNames[2:length(dfNames)], reorderData, variantList = variantOrder))

## Concatenate refCount/altCount columns according to donor/condition order 
# Find and separate dfNames based on condition order
condition1_dfNames <- dfNames[str_detect(dfNames, conditions[1])]
condition2_dfNames <- dfNames[str_detect(dfNames, conditions[2])]

# Put each condition-separated dfNames into order based on donor order (this will also automatically order by bio rep)
donor_condition1_order <- unlist(lapply(donors, donorMatch, condition1_dfNames))
donors_condition1_dfNames <- condition1_dfNames[donor_condition1_order]
donor_condition2_order <- unlist(lapply(donors, donorMatch, condition2_dfNames))
donors_condition2_dfNames <- condition2_dfNames[donor_condition2_order]

## Assemble and cbind one dataframe, alternating between grabbing refCount/altCount columns from dfs in conditions
alleleCounts <- data.frame(matrix(nrow = length(variantOrder)))

# Make corresponding sample table for later use in DESeq2 (rows correspond to columns in alleleCounts) (can add more info to this later) 
colData <- data.frame(matrix(ncol = 4))
colnames(colData) <- c("Donor", "Treatment", "BioRep", "Allele")

for (i in 1:max(length(donors_condition1_dfNames), length(donors_condition2_dfNames))){


    # Grab condition1 data and set new column names based on donor/condition/biorep
    if (!is.na(donors_condition1_dfNames[i])){
        cond1_data <- get(donors_condition1_dfNames[i])
        cond1_data <- cond1_data[,c("refCount", "altCount")]
        donor <- unlist(str_split(donors_condition1_dfNames[i], "_"))[2]
        cond1_biorep <- unlist(str_split(donors_condition1_dfNames[i], "_"))[5]
        colnames(cond1_data) <- c(paste(donor, conditions[1], cond1_biorep, "ref", sep = "_"), paste(donor, conditions[1], cond1_biorep, "alt", sep = "_"))
        alleleCounts <- cbind(alleleCounts, cond1_data)
        colData <- rbind(colData, c(donor, conditions[1], cond1_biorep, "ref"),
                                c(donor, conditions[1], cond1_biorep, "alt"))
    }
    
    if (!is.na(donors_condition2_dfNames[i])){
        # Grab condition2 data and set new column names based on donor/condition/biorep
        cond2_data <- get(donors_condition2_dfNames[i])
        cond2_data <- cond2_data[,c("refCount", "altCount")]
        donor <- unlist(str_split(donors_condition2_dfNames[i], "_"))[2]
        cond2_biorep <- unlist(str_split(donors_condition2_dfNames[i], "_"))[5]
        colnames(cond2_data) <- c(paste(donor, conditions[2], cond2_biorep, "ref", sep = "_"), paste(donor, conditions[2], cond2_biorep, "alt", sep = "_"))
        alleleCounts <- cbind(alleleCounts, cond2_data)
        colData <- rbind(colData, c(donor, conditions[2], cond2_biorep, "ref"),
                            c(donor, conditions[2], cond2_biorep, "alt"))
    }
        
} 
# Remove first column/row of NA's
alleleCounts <- alleleCounts[,-1]
colData <- colData[-1,]

# Set rownames in alleleCounts as snps
rownames(alleleCounts) <- variantOrder

## Write allele counts matrix and sample data to output
write.table(alleleCounts, "output/AI/alleleCountsMatrix.txt", quote = FALSE, row.names = TRUE, col.names = TRUE)
write.table(colData, "output/AI/colData.txt", quote = FALSE, row.names = FALSE, col.names = TRUE)
