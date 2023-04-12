library(data.table)
library(DESeq2)
library(googlesheets4)
library(dplyr)

# colData ######################################################################

colData <- fread("../processing/output/AI/colData.csv", data.table = FALSE)


# donorInfo <- read_sheet(ss = "https://docs.google.com/spreadsheets/d/1JwLw9D6rMqhHC9BPrZebAN40Wojo-CqbMdiIPXAzkLo/edit#gid=1699779981",
#                         sheet = "Donors") %>% dplyr::rename(donor = Donor)
# donorInfo$donor <- unlist(donorInfo$donor)

# colData <- colData %>% left_join(donorInfo[,c("donor", "Sex", "Age")])

# Allele counts matrix #########################################################

alleleCountsMatrix <- fread("../processing/output/AI/alleleCountsMatrix.csv", data.table = FALSE,
                            colClasses = c("character", rep("numeric", nrow(colData))),
                            skip = 4)

# Set rownames as first column (variant IDs)
rownames(alleleCountsMatrix) <- alleleCountsMatrix$V1

# Remove variantID first column
alleleCountsMatrix <- alleleCountsMatrix[,-1]

# Coerce to matrix
alleleCountsMatrix <- as.matrix(alleleCountsMatrix)

# Weight matrix ################################################################

weightMatrix <- fread("../processing/output/AI/weightsMatrix.csv", skip = 4, data.table = FALSE,
                      colClasses = c("character", rep("numeric", nrow(colData))))
rownames(weightMatrix) <- weightMatrix$V1
weightMatrix <- weightMatrix[,-1]
weightMatrix <- as.matrix(weightMatrix)




# DESeq ########################################################################
design <- ~0 + condition:donor + condition:allele

# Add datasets, weights, and size factors
dds <- DESeqDataSetFromMatrix(alleleCountsMatrix, colData, design)
assays(dds)[["weights"]] <- weightMatrix
# No size factor normalization, set all to 1
sizeFactors(dds) <- rep(1, ncol(alleleCountsMatrix))
# Relevel with reference allele set first
dds$allele <- relevel(dds$allele, ref = "ref")
save(dds, file = "data/dds.rda")
allelic_imbalance_analysis <- DESeq(dds, minReplicatesForReplace = Inf, fitType = "local")
save(allelic_imbalance_analysis, file = paste0("data/", Sys.Date(), "_AIanalysis.rda"))

resCTL = results(allelic_imbalance_analysis,name = "conditionCTL.allelealt")
notNA_resCTL = resCTL[!is.na(resCTL$log2FoldChange),]
save(notNA_resCTL, file = paste0("data/", Sys.Date(), "_AIresCTL.rda"))
resFNF = results(allelic_imbalance_analysis,name = "conditionFNF.allelealt")
notNA_resFNF = resFNF[!is.na(resFNF$log2FoldChange),]
save(notNA_resFNF, file = paste0("data/", Sys.Date(), "_AIresFNF.rda"))
