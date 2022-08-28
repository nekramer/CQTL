library(data.table)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(dplyr)
source('../AIanalysis/scripts/utils.R')

args = commandArgs(trailingOnly = TRUE)
# args[1] is the ASEP data
# args[2] is the chromosome
# args[3] is the line number

getGenes <- function(countData_chunk, chrom, chunkNum, txdb = TxDb.Hsapiens.UCSC.hg38.knownGene, orgdb = org.Hs.eg.db){
    # Parse 'snp' column and convert to GRanges
    ranges <- varID_to_GRanges(countData_chunk$snp)

    # Get genes from ranges
    genes <- GRanges_to_Genes(ranges, txdb = txdb, orgdb = orgdb, singleStrandOnly = FALSE)
    
    # Join back to countData
    countData_ranges <- cbind(countData_chunk, as.data.frame(ranges)[,c("seqnames", "start")])
    colnames(countData_ranges) <- c("id", "group", "snp", "ref", "total", "chr", "pos")
    ASEP_genes <- left_join(countData_ranges, genes)

    # Only keep necessary columns and rename gene_symbol to gene
    ASEP_genes_final <- ASEP_genes[,c("id", "group", "snp", "ref", "total", "gene_symbol")]
    colnames(ASEP_genes_final)[6] <- c("gene")

    # Keep rest of gene info in separate dataframe
    gene_info <- ASEP_genes[,c("gene_symbol", "gene_start", "gene_end", "gene_strand")]

    # Write to files
    fwrite(ASEP_genes_final, file = paste0("output/AI/chr", chrom, "_ASEPfinal_chunk", chunkNum, ".txt"), col.names = FALSE)
    fwrite(gene_info, file = paste0("output/AI/chr", chrom, "_ASEPfinal_geneInfo_chunk", chunkNum,".txt"), col.names = FALSE)
}

# Get header columns
header <- colnames(fread(args[1], nrows = 0))

# Calculate number of chunks based on line number (-1 for header)
nChunks <- ceiling((as.numeric(args[3]) - 1)/1000)

skip <- 0
for (chunk in 1:nChunks){
  dataChunk <- fread(args[1], nrows = 100, skip = skip, data.table = FALSE)
  
  # Assign colnames to chunk
  colnames(dataChunk) <- header
  
  # Get genes and write file for that chunk
  getGenes(dataChunk, chrom = args[2], chunkNum = chunk)
  
  # Increment skip to get next chunk
  skip <- skip + 100
}

# Write files for headers
fwrite(data.frame(matrix(header, nrow = 1)), file = paste0("output/AI/chr", args[2], "_ASEPheader.txt"), col.names = FALSE)
fwrite(data.frame(matrix(c("gene_symbol", "gene_start", "gene_end", "gene_strand"), nrow = 1)), file = paste0("output/AI/chr", args[2], "_geneheader.txt"), col.names = FALSE)