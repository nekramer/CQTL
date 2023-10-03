library(DESeq2)
library(tidyverse)
library(plyranges)

# Load gse object
load("data/2023-10-03_gse.rda")

# Read in donorSamplesheet for additional donor info
donorSamplesheet <- read_csv("data/donorSamplesheet.csv") |> 
  mutate(Race = replace_na(Race, "Unknown"))

# Join gse colData with donorSamplesheet
colData(gse) <- as(left_join(as.data.frame(colData(gse)),
                             donorSamplesheet[,c("Donor", "Sex", "Age", "Race")],
                             by = "Donor"), "DataFrame")
# Convert colData to factors
colData(gse)[] <- lapply(colData(gse), factor)

# Build DESeq object
dds <- DESeqDataSet(gse, design = ~Donor + Condition)
colnames(dds) <- colData(gse)[,"names"]

# Filter out lowly expressed genes
keep <- rowSums(counts(dds) >= 10) >= ceiling(nrow(colData(gse))*0.10)
dds <- dds[keep,]

# Fit model
dds <- DESeq(dds)

## Save dds
save(dds, file = "data/condition_de/differential_expression_dds.rda")

# Shrink l2fc
de_genes_shrink <- lfcShrink(dds,
                             coef = "Condition_FNF_vs_CTL", format = "GRanges") |>
  plyranges::names_to_column("gene_id")

# Join results with gene info
de_genes_shrink <-
  inner_join(x = as.data.frame(de_genes_shrink),
             y = as.data.frame(rowData(gse)) %>%
               dplyr::select(c("gene_id", "symbol", "tx_ids")),
             by = "gene_id") %>%
  makeGRangesFromDataFrame(keep.extra.columns = TRUE) %>%
  keepStandardChromosomes(pruning.mode = "coarse") %>%
  as.data.frame()

## Save l2fc-shrunken results
write_csv(de_genes_shrink, file = "data/condition_de/de_genes_results.csv")

# Get significant genes
sig_deGenes_pval01_l2fc2 <- de_genes_shrink |> 
  filter(padj < 0.01 & abs(log2FoldChange) > 2) |> 
  write_csv("data/condition_de/sig_deGenes_pval01_l2fc2.csv")

# Split into upregulated and downregulated
upsig_deGenes_pval01_l2fc2 <- sig_deGenes_pval01_l2fc2 |> 
  filter(log2FoldChange > 0) |> 
  write_csv("data/condition_de/upsig_deGenes_pval01_l2fc2.csv")
downsig_deGenes_pval01_l2fc2 <- sig_deGenes_pval01_l2fc2 |> 
  filter(log2FoldChange < 0) |> 
  write_csv("data/condition_de/downsig_deGenes_pval01_l2fc2.csv")


# Run Homer for GO terms, KEGG pathways, and TF binding motifs ------------

# Upregulated genes
system("scripts/run_homer.sh data/condition_de/upsig_deGenes_pval01_l2fc2.csv data/homer/homer_upsig_deGenes_pval01_l2fc2")

# Downregulated genes
system("scripts/run_homer.sh data/condition_de/downsig_deGenes_pval01_l2fc2.csv data/homer/homer_downsig_deGenes_pval01_l2fc2")
