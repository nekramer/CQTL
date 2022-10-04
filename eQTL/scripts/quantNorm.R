#!/usr/bin/R
library(tximeta)
library(readr)
library(dplyr)
library(tibble)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(DESeq2)
library(edgeR)

# FUNCTIONS -------------------------------------------------------------------
# Function to inverse normalize a row of gene counts
inverseNormGene <- function(geneRow){
  normValues <- qnorm((rank(as.numeric(geneRow),
                            na.last = "keep") - 0.5)/sum(!is.na(as.numeric(geneRow))))
  return(normValues)
}
# READ IN ---------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
samplesheet <- args[1]
sampleQuants <- args[2:length(args)]
samples <- basename(dirname(sampleQuants))

coldata <- read_csv(samplesheet) %>%
  # Get distinct samples
  distinct(Sample, .keep_all = TRUE) %>%
  # Put in order of input quants
  arrange(match(Sample, samples)) %>%
  # add files
  mutate(files = sampleQuants) %>%
  # names column
  mutate(names = Sample) %>%
  # Condition groups
  mutate(Condition = as.factor(Condition)) %>% 
  mutate(group = as.numeric(Condition))

# Import salmon transcript quantification-------------------------------------
se <- tximeta(coldata)

# Convert to gene-level scaled transcripts -------------------------------------
gse <- summarizeToGene(se)

# Filter out lowly expressed genes ---------------------------------------------
# 5 reads in at least 25% of samples?
# 10 reads in at least 5% of samples?
keep <- filterByExpr(gse, min.count = 10, min.prop = 0.05)
gse_filtered <- gse[keep,]

# TMM normalization ---------------------------------------------------------
gse_quant <- calcNormFactors(gse_filtered, method = "TMM")

# Grab CPM counts ---------------------------------------------------------
CQTL_CPMadjTMM <- as.data.frame(cpm(gse_quant))

# Inverse normal transformation -----------------------------------------------
# Split based on condition
CTL_CPMadjTMM <- dplyr::select(CQTL_CPMadjTMM, contains("CTL"))
FNF_CPMadjTMM <- dplyr::select(CQTL_CPMadjTMM, contains("FNF"))

# Inverse normalize across genes in each condition
CTL_CPMadjTMM_invNorm <- as.data.frame(t(apply(CTL_CPMadjTMM, 1, inverseNormGene))) %>%
  rownames_to_column("gene_id")
colnames(CTL_CPMadjTMM_invNorm) <- c("gene_id", colnames(CTL_CPMadjTMM))
FNF_CPMadjTMM_invNorm <- as.data.frame(t(apply(FNF_CPMadjTMM, 1, inverseNormGene))) %>%
  rownames_to_column("gene_id")
colnames(FNF_CPMadjTMM_invNorm) <- c("gene_id", colnames(FNF_CPMadjTMM))

# Also inverse normalize across genes in combined conditions
CPMadjTMM_invNorm <- as.data.frame(t(apply(CQTL_CPMadjTMM, 1, inverseNormGene))) %>%
  rownames_to_column("gene_id")
colnames(CPMadjTMM_invNorm) <- c("gene_id", colnames(CQTL_CPMadjTMM))

# Join with gene info, add length column, and reorder into BED format --------
gene_info <- as.data.frame(rowRanges(gse_filtered)) %>% 
  dplyr::select(seqnames, start, end, strand, gene_id, gene_name) %>%
  mutate(seqnames = paste0("chr", seqnames))

CTL_CPMadjTMM_invNorm <- CTL_CPMadjTMM_invNorm %>% left_join(gene_info) %>%
  relocate(seqnames, start, end, gene_id, gene_name, strand) %>%
  rename("seqnames" = "#chr")
FNF_CPMadjTMM_invNorm <- FNF_CPMadjTMM_invNorm %>% left_join(gene_info) %>%
  relocate(seqnames, start, end, gene_id, gene_name, strand) %>%
  rename("seqnames" = "#chr")
CPMadjTMM_invNorm <- CPMadjTMM_invNorm %>% left_join(gene_info) %>%
  relocate(seqnames, start, end, gene_id, gene_name, strand) %>%
  rename("seqnames" = "#chr")

write_delim(CTL_CPMadjTMM_invNorm, file = "output/normquant/CTL_CPMadjTMM_invNorm.bed", delim = "\t")
write_delim(FNF_CPMadjTMM_invNorm, file = "output/normquant/FNF_CPMadjTMM_invNorm.bed", delim = "\t")
write_delim(CPMadjTMM_invNorm, file = "output/normquant/ALL_CPMadjTMM_invNorm.bed", delim = "\t")

