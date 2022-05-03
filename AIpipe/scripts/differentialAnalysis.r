#!/usr/bin/bash Rscript
library(data.table)
library(DESeq2)
library(stringr)

args = commandArgs(trailingOnly = TRUE)
# args[1]: alleleCounts matrix
# args[2]: colData for alleleCounts matrix
# args[3]: weight matrix

# Read in filtered alleleCounts matrix
alleleCountsMatrix <- read.table(args[1])

# Convert rownames to separate dataframe of chrom and position
positions <- data.frame(do.call(rbind, str_split(rownames(alleleCountsMatrix), ":")))
positions <- positions[,c(1,2)]
colnames(positions) <- c("chr", "pos")
write.table(positions, file = "output/AI/AI_variant_positions.txt", quote = FALSE, row.names = FALSE, col.names = TRUE)

# Read in sampleData 
colData <- fread(args[2], data.table = FALSE)

# Read in weight matrix and set row names
weightMatrix <- fread(args[3], data.table = FALSE)
rownames(weightMatrix) <- weightMatrix$V1
weightMatrix <- as.matrix(weightMatrix[,-1])

# Initialize design for DESeq
design <- ~0 + Treatment:Donor + Treatment:Allele

# Add datasets, weights, and size factors
dds <- DESeqDataSetFromMatrix(alleleCountsMatrix, colData, design)
assays(dds)[["weights"]] <- weightMatrix
# No size factor normalization, set all to 1
sizeFactors(dds) <- rep(1, ncol(alleleCountsMatrix))
# Relevel with reference allele set first
dds$Allele <- relevel(dds$Allele, ref = "ref")

allelic_imbalance_analysis <- DESeq(dds)

# Save allelic imbalance object
save(allelic_imbalance_analysis, file = "output/AI/differentialAllelicImbalance.rda")




# pullgenedata <- function(gene,dds,weightMatrix,sample_Data)
# {
#   # pull the data
#   a= plotCounts(gene=gene,dds = dds,intgroup  = c("condition","allele","donor"),returnData=TRUE)
  
#   # add the donor just to make sure it agrees with original order
#   a$test = sample_Data$donor
  
#   # pull the weights
#   a$het = weightMatrix[gene,] > .5
  
#   # filter for hets
#   a = a[a$het==TRUE,]
  
#   # remove the sample with 2 reps
#   if ("AM7205" %in% a$donor)
#   {
#     a = a[-which(a$donor == "AM7205"),]
#   }
  
#   # make a new cleaned up df
#   b = data.frame(donor=a$donor,condition=a$condition,allele=a$allele,count=a$count)
  
#   # convert to wide format
#   data_wide <- pivot_wider(data=b, names_from = c(condition,allele),values_from = count)
  
#   # get rid of NAs probably not hets
#   data_wide = data_wide[complete.cases(data_wide),]
  
#   # organize columns
#   if (nrow(data_wide) > 0 & ncol(data_wide) ==5)
#   {
#     data_wide = data_wide[,c("donor","CTL_ref","CTL_alt","FNF_ref","FNF_alt")]
#   }
#   return (data_wide)
# }

# rsidsToGenes <- function(rsids, snpDb, TxDb, orgDb){
#   snpDb <- eval(parse(text = snpDb))
#   TxDb <- eval(parse(text = TxDb))
#   orgDb <- eval(parse(text = orgDb))
#   ## Grab snp locations based on rsid
#   snp_locations <- snpsById(snpDb, ids = rsids, ifnotfound = "drop")
#   seqlevelsStyle(snp_locations) <- "UCSC"
#   ## Grab overall set of genes that overlap with the snp locations
#   genes <- suppressMessages(genes(TxDb))
#   genes_snp_subset <- subsetByOverlaps(genes, snp_locations)
#   ## Find matching keys of overlaps
#   gene_snp_overlaps <- findOverlaps(genes_snp_subset, snp_locations)
#   ## Get gene names fromo OrgDb
#   gene_overlaps <- genes_snp_subset[from(gene_snp_overlaps),]
#   gene_symbols <- suppressMessages(AnnotationDbi::select(orgDb, keys = gene_overlaps$gene_id, columns = "SYMBOL", keytype = "ENTREZID"))
#   gene_overlaps_symbols <- merge(gene_overlaps, gene_symbols, by.x = "gene_id", by.y = "ENTREZID")
#   ## Combine matched gene names with matched snps by matching geneids
#   snp_overlaps <- cbind(gene_overlaps$gene_id, as.data.frame(snp_locations[to(gene_snp_overlaps),]))
#   colnames(snp_overlaps)[1] <- "gene_id"
#   snp_genes <- dplyr::left_join(snp_overlaps, gene_overlaps_symbols, by = "gene_id")
#   snp_genes <- snp_genes[,c(1:3, 5, 8:9, 11:12)]
#   snp_genes <- unique(snp_genes)
#   colnames(snp_genes) <- c("gene_ENTREZID", "chrom", "pos", "rsid", "gene_start", "gene_end", "gene_strand", "gene_SYMBOL")
#   ## slightly reorder columns for side-by-side viewing of rsids and gene symbols
#   snp_genes <- snp_genes[c(2:3, 5:7, 1, 8, 4)]
#   return(snp_genes)
# }