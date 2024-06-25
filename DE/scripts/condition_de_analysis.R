library(DESeq2)
library(tidyverse)
library(plyranges)
library(googledrive)
library(googlesheets4)

# Load gse object
load("data/2023-10-03_gse.rda")
# Remove geno-contaminated donor
gse <- gse[, gse$Donor != "AM7352"]
# Read in donorSamplesheet for additional donor info
donorSamplesheet <- read_csv("data/donorSamplesheet.csv") |> 
  mutate(Race = replace_na(Race, "Unknown")) |> 
  dplyr::select(Donor, Sex, Age) |> 
  # Read in and join ancestries determined through genotyping pca
  left_join(read_csv("data/CQTL_COA_01_GDA8_COA2_01_COA3_01_GDA8_COA4_COA5_COA6_COA7_predictedAncestry.csv") |> 
              separate_wider_delim(cols = "Donor", delim = "_", 
                                   names = c(NA, "Donor", NA), too_many = "drop"), by = "Donor")

# Join gse colData with donorSamplesheet
colData(gse) <- as(left_join(as.data.frame(colData(gse)),
                             donorSamplesheet,
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


## Reformat and write to supplementary table
sup_table_de <- de_genes_shrink |> 
  dplyr::select(symbol, gene_id, padj, log2FoldChange) |> 
  mutate(`FNF response` = case_when(padj < 0.01 &
                                      abs(log2FoldChange) > 2 &
                                      log2FoldChange < 0 ~ "---",
                                    padj < 0.01 &
                                      abs(log2FoldChange) > 2 &
                                      log2FoldChange > 0 ~ "+++", 
                                    padj < 0.05 &
                                      abs(log2FoldChange) > 1 &
                                      log2FoldChange < 0 ~ "-",
                                    padj < 0.05 & 
                                      abs(log2FoldChange) > 1 &
                                      log2FoldChange > 0 ~ "+")) |> 
  filter(!is.na(`FNF response`)) |> 
  arrange(symbol)

write_csv(sup_table_de, file = "tables/SupTable1.csv")

# Write to google drive
ss <- gs4_create(name = "SupTable1")
write_sheet(sup_table_de,
            ss, sheet = "Sheet1")
drive_mv(file = "SupTable1", path = as_dribble("CQTL paper/Figures and Tables"))

 # Significant genes at various thresholds ---------------------------------

sig_deGenes_pval05_l2fc1 <- de_genes_shrink |> 
  filter(padj < 0.05 & abs(log2FoldChange) > 1) |> 
  write_csv("data/condition_de/sig_deGenes_pval05_l2fc1.csv")

sig_deGenes_pval01_l2fc1 <- de_genes_shrink |> 
  filter(padj < 0.01 & abs(log2FoldChange) > 1) |> 
  write_csv("data/condition_de/sig_deGenes_pval01_l2fc1.csv")

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
