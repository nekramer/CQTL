library(data.table)
library(DESeq2)

# Allele counts matrix #########################################################

alleleCountsMatrix <- fread("alleleCountsMatrix.csv", data.table = FALSE,
                            colClasses = c("character", rep("numeric", 224)),
                            skip = 4)

# Set rownames as first column (variant IDs)
rownames(alleleCountsMatrix) <- alleleCountsMatrix$V1

# Remove variantID first column
alleleCountsMatrix <- alleleCountsMatrix[,-1]

# Coerce to matrix
alleleCountsMatrix <- as.matrix(alleleCountsMatrix)

# Weight matrix ################################################################

weightMatrix <- fread("weightsMatrix.csv", skip = 4, data.table = FALSE,
                      colClasses = c("character", rep("numeric", 224)))
rownames(weightMatrix) <- weightMatrix$V1
weightMatrix <- weightMatrix[,-1]
weightMatrix <- as.matrix(weightMatrix)

# colData ######################################################################

colData <- fread("colData.csv", data.table = FALSE)

# DESeq ########################################################################
design <- ~0 + condition:donor + condition:allele

# Add datasets, weights, and size factors
dds <- DESeqDataSetFromMatrix(alleleCountsMatrix, colData, design)
assays(dds)[["weights"]] <- weightMatrix
# No size factor normalization, set all to 1
sizeFactors(dds) <- rep(1, ncol(alleleCountsMatrix))
# Relevel with reference allele set first
dds$allele <- relevel(dds$allele, ref = "ref")

allelic_imbalance_analysis <- DESeq(dds)
