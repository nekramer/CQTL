#!/usr/bin/env Rscript

library(data.table)
library(stringr)
library(dplyr)
library(purrr)
library(pryr)
library(bigmemory)
args = commandArgs(trailingOnly = TRUE)
# args[1]: donors
# args[2]: conditions
# args[3]...args[n]: alleleCount files from each sample

## Functions
# Function to read in and assign data to dataframe
readData <- function(file){
    name <- gsub("_alleleCounts.csv", "", basename(file))
    data <- fread(file, data.table = FALSE)
    data_subset <- data[,c("variantID", "refCount", "altCount")]

    # Rename refCount and altCount columns based on sample name
    colnames(data_subset)[2:3] <- c(paste0(name, "_ref"), paste0(name, "_alt"))
    return(data_subset)
}

# Function to parse BioRep slot in sample names
getBioRep <- function(file){
    name <- gsub("_alleleCounts.csv", "", basename(file))
    biorep <- unlist(str_split(name, "_"))[6]
    return(biorep)
}

## Read in all data
allData <- lapply(args[3:length(args)], readData)
print("read data")

## Full join all the data together and set first column as rownames
alleleCountsMatrix <- as.matrix(purrr::reduce(allData, dplyr::full_join, by = "variantID"))
print("joined")
#print(paste0("Object size: ", object_size(alleleCountsMatrix)))

# options(bigmemory.allow.dimnames = TRUE)
rownames(alleleCountsMatrix) <- alleleCountsMatrix[,"variantID"]
alleleCountsMatrix <- alleleCountsMatrix[,-1]
print("renamed row")
print(paste0("Mem used1: ", mem_used()))

## Initialize weight matrix with same size with all 1s
# weightsMatrix <- matrix(data = 1, nrow = nrow(alleleCountsMatrix), ncol = ncol(alleleCountsMatrix))
# rownames(weightsMatrix) <- rownames(alleleCountsMatrix)
# colnames(weightsMatrix) <- colnames(alleleCountsMatrix)
weightsMatrix <- big.matrix(init = 1, nrow = nrow(alleleCountsMatrix), ncol = ncol(alleleCountsMatrix),
                            dimnames = list(rownames(alleleCountsMatrix), colnames(alleleCountsMatrix)))
print("initialized weight matrix")
print(paste0("Mem used2: ", mem_used()))
indx <- which(is.na(alleleCountsMatrix), arr.ind = TRUE)
weightMatrix[indx] <- 0
print(paste0("Mem used3: ", mem_used()))

## Make slots of weightsMatrix where alleleCountsMatrix is NA a 0
# weightsMatrix[is.na(alleleCountsMatrix)] <- 0

## Now convert NAs to 0s in alleleCountsMatrix
alleleCountsMatrix[is.na(alleleCountsMatrix)] <- 0

## Create corresponding sample data 
donors <- unlist(str_split(args[1], ","))
conditions <- unlist(str_split(args[2], ","))
alleles <- c("ref", "alt")

colData <- expand.grid(donors, conditions, alleles) %>%
    arrange(factor(Var1, levels = donors),
            factor(Var2, levels = conditions),
            factor(Var3, levels = alleles)) %>%
    rename(Donor = Var1, Treatment = Var2, Allele = Var3)

bioreps <- unlist(lapply(args[3:length(args)], getBioRep))
colData$BioRep <- bioreps

## Write alleleCountsMatrix, weightsMatrix, and colData to output files
#write.table(alleleCountsMatrix, "output/AI/alleleCountsMatrix.txt", quote = FALSE, row.names = TRUE, col.names = TRUE)
#write.table(weightsMatrix, "output/AI/weightsMatrix.txt", quote = FALSE, row.names = TRUE, col.names = TRUE)
#write.table(colData, "output/AI/colData.txt", quote = FALSE, row.names = FALSE, col.names = TRUE)