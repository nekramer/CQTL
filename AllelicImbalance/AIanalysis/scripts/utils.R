# Function to parse variant IDs in the format chrX:pos:ref:alt and
# convert to GRanges/dataframe
library(stringr)
library(purrr)
library(plyranges)
varID_to_GRanges <- function(varIDs, GRanges = TRUE){
  
  chr <- varIDs %>% str_split(":") %>% map(1) %>% unlist()
  pos <- varIDs %>% str_split(":") %>% map(2) %>% unlist()
  
  ranges <- GRanges(seqnames = chr, 
                    ranges = IRanges(start = as.numeric(pos), 
                                     end = as.numeric(pos)))
  if (GRanges == FALSE){
    ranges <- as.data.frame(ranges)
  }
  return(ranges) 
}

GRanges_to_Genes <- function(ranges, txdb, orgdb, singleStrandOnly = TRUE){
  genes <- subsetByOverlaps(suppressMessages(genes(txdb, 
                                  single.strand.genes.only = singleStrandOnly)), 
                            ranges)
  
  if (singleStrandOnly == FALSE){
    
    genes <- bind_ranges(as.list(genes), .id = "gene_id")
    genes$gene_id <- as.character(genes$gene_id)
  }
  
  # Find matching keys of overlaps and grab those genes
  overlaps <- findOverlaps(genes, ranges)
  
  gene_overlaps <- genes[from(overlaps)]
  snp_overlaps <- cbind(gene_overlaps$gene_id, as.data.frame(ranges[to(overlaps),]))
  colnames(snp_overlaps)[1] <- "gene_id"

  # Get gene names from orgDb
  ## Get gene id type from txdb metadata
  txdb_metadata <- metadata(txdb)
  geneid <- tolower(txdb_metadata[txdb_metadata$name == "Type of Gene ID",2])
  
  if (length(grep("entrez", geneid)) == 1){
    gene_keytype <- "ENTREZID"
  } else if (length(grep("ensembl", geneid)) == 1){
    gene_keytype <- "ENSEMBL"
    # Remove decimal of Ensembl gene id
    gene_overlaps$gene_id <- gsub("\\..*", "", gene_overlaps$gene_id)
    snp_overlaps$gene_id <- gsub("\\..*", "", snp_overlaps$gene_id)
  }
  
  gene_symbols <- suppressMessages(AnnotationDbi::select(orgdb, keys = gene_overlaps$gene_id, columns = "SYMBOL", keytype = gene_keytype))
  # Recombine with original gene information
  gene_overlaps_symbols <- merge(gene_overlaps, gene_symbols, by.x = "gene_id", by.y = gene_keytype)
  # Combine matched gene names with matched snps by matching geneids
  snp_genes <- left_join(snp_overlaps, gene_overlaps_symbols, by = "gene_id")
  # Get unique combinations
  snp_genes <- unique(snp_genes)
  
  # Filter and reorganize 
  snp_genes <- snp_genes[,c(2, 3, 12, 8, 9, 11, 1)]
  colnames(snp_genes) <- c("chr", "pos", "gene_symbol", "gene_start", "gene_end", "gene_strand", gene_keytype)
  
  return(snp_genes)
}
