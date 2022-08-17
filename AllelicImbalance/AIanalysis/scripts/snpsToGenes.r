#!/usr/bin/env Rscript

library(optparse)


option_list <- list(
    make_option(c("--file"), type = "character", default = NULL,
                help = "Input file of snps with location and rsid (chr column, pos column, rsid column)."),
    make_option(c("--txdb"), type = "character", default = "TxDb.Hsapiens.UCSC.hg38.knownGene",
                help = "Name of TxDb to use for gene transcript information or a path to a custom sqlite TxDb. (default %default)"),
    make_option(c("--orgdb"), type = "character", default = "org.Hs.eg.db",
                help = "Name of orgDb to use for gene symbol information. (default %default)"),
    make_option(c("--output"), type = "character", default = NULL,
                help = "Name of desired output file.")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(AnnotationDbi))
suppressPackageStartupMessages(library(tools))

# Read in variant file locations
variants <- fread(opt$file, data.table = FALSE)

# Make sure this file has 3 columns, each named chr, pos, rsid
if (ncol(variants) != 3){
    stop("Input snp file must have 3-columns specifying the chromosome, position, and rsid.")
}
if (!"chr" %in% colnames(variants) | !"pos" %in% colnames(variants) | !"rsid" %in% colnames(variants)){
    stop("Input snp file must have columns labeled `chr`, `pos`, and `rsid`.")
}


# Convert positions to GRanges
variant_ranges <- GRanges(seqnames = variants$chr, ranges = IRanges(start = variants$pos, end = variants$pos), rsid = variants$rsid)

# Load TxDb (check if created or included library)
if (file_ext(opt$txdb) == "sqlite"){
    txdb <- loadDb(opt$txdb)
} else {
    suppressPackageStartupMessages(library(opt$txdb, character.only = TRUE))
    txdb <- eval(parse(text = opt$txdb))
}

# Load orgDb
suppressPackageStartupMessages(library(opt$orgdb, character.only = TRUE))
orgdb <- eval(parse(text = opt$orgdb))

# Grab set of genes that overlap with variant locations
genes <- subsetByOverlaps(genes(txdb), variant_ranges)

# Find matching keys of overlaps and grab those genes
overlaps <- findOverlaps(genes, variant_ranges)
gene_overlaps <- genes[from(overlaps)]
snp_overlaps <- cbind(gene_overlaps$gene_id, as.data.frame(variant_ranges[to(overlaps),]))
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
snp_genes <- snp_genes[,c(7, 2, 3, 13, 9, 10, 12, 1)]
colnames(snp_genes) <- c("rsid", "chr", "pos", "gene_symbol", "gene_start", "gene_end", "gene_strand", gene_keytype)

# Sort by chromosome
snp_genes$chr <- as.factor(snp_genes$chr)
snp_genes <- snp_genes %>% arrange(chr)

# Write to output
write.table(snp_genes, file = opt$output, quote = FALSE, row.names = FALSE, col.names = TRUE)
