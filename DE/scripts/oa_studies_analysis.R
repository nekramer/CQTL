library(GEOquery)
library(DESeq2)
library(readxl)


## GSE114007
getGEOSuppFiles("GSE114007", baseDir = "data/GSE114007/", makeDirectory = FALSE)

# Read in raw counts for normal and OA from GSE114007_raw_counts.xlsx
GSE114007_normal_raw_counts <- read_excel("data/GSE114007/GSE114007_raw_counts.xlsx",
                                          sheet = "Normal")
GSE114007_oa_raw_counts <- read_excel("data/GSE114007/GSE114007_raw_counts.xlsx",
                                          sheet = "OA")

# Join all raw counts
GSE114007_raw_counts <- full_join(GSE114007_normal_raw_counts, 
                                  GSE114007_oa_raw_counts, by = "symbol")

# Make colData and count matrix for DESeq
GSE114007_colData <- tibble(Sample = colnames(GSE114007_raw_counts)[-1]) |> 
  # Batch
  # mutate(Batch = case_when(Sample %in% c("Normal_Cart_2_2", "Normal_Cart_3_3",
  #                                        "Normal_Cart_4_4", "Normal_Cart_5_5",
  #                                        "Normal_Cart_6_6", "Normal_Cart_7_3",
  #                                        "Normal_Cart_9_7", "Normal_Cart_10_8",
  #                                        "OA_Cart_1_7", "OA_Cart_2_8", "OA_Cart_3_9",
  #                                        "OA_Cart_4_10", "OA_Cart_5_5", "OA_Cart_6_1",
  #                                        "OA_Cart_7_2", "OA_Cart_8_5", "OA_Cart_9_6",
  #                                        "OA_Cart_10_9") ~ "Batch1",
  #                          Sample %in% c("normal_01", "normal_02", "normal_03", "normal_04",
  #                                        "normal_05", "normal_06", "normal_07", "normal_08",
  #                                        "normal_09", "normal_10", "OA_01", "OA_02", 
  #                                        "OA_03", "OA_04", "OA_05", "OA_06", "OA_07",
  #                                        "OA_08", "OA_09", "OA_10") ~ "Batch2")) |> 
  # Reformat strings with "Cart"
  mutate(Sample = gsub("_Cart", "", Sample)) |> 
  # Collapse numbers when split by _
  mutate(Sample = gsub("_([0-9]+)_([0-9]+)", "_\\1\\2", Sample)) |> 
  # Split into Condition and a dummy "donor number"
  separate_wider_delim(cols = "Sample",
                       delim = "_",
                       names = c("Condition", "Donor")) |> 
  # Convert Condtions to all lowercase for consistency
  mutate(Condition = tolower(Condition)) |> 
  # Condition factors
  mutate(Condition = factor(Condition, levels = c("normal", "oa")))

GSE114007_countMatrix <- GSE114007_raw_counts |> 
  column_to_rownames(var = "symbol") |> 
  as.matrix() 

# Make DESeqDataSet
GSE114007_dds <- DESeqDataSetFromMatrix(GSE114007_countMatrix, 
                                        GSE114007_colData, design = ~Condition) 

# Filter out lowly expressed genes
keep <- rowSums(counts(GSE114007_dds) >= 10) >= 4
GSE114007_dds <- GSE114007_dds[keep,]

# Fit model
GSE114007_dds <- DESeq(GSE114007_dds)
save(GSE114007_dds, file = "data/GSE114007/GSE114007_dds.rda")

# Shrink l2fc
GSE114007_shrink <- lfcShrink(GSE114007_dds,
                             coef = "Condition_oa_vs_normal") |>
  as_tibble(rownames = "symbol")

# All gene results
write_csv(GSE114007_shrink, file = "data/GSE114007/GSE114007_deseq_res.csv")

# Get significant genes; this study used padj < 0.05 and abs(log2FoldChange) > 1 significance cutoffs
GSE114007_pval05_l2fc1 <- GSE114007_shrink |> 
  filter(padj < 0.05 & abs(log2FoldChange) > 1) |> 
  write_csv("data/GSE114007/GSE114007_pval05_l2fc1.csv")

## GSE168505
getGEOSuppFiles("GSE168505", baseDir = "data/GSE168505/", makeDirectory = FALSE)

# Read in count matrix
GSE168505_raw_counts <- read_tsv("data/GSE168505/GSE168505_osteoarthritis_matrix.tsv.gz")

# Make colData and count matrix for DESeq
GSE168505_colData <- tibble(Sample = colnames(GSE168505_raw_counts)[-1]) |> 
  separate_wider_delim(cols = "Sample", delim = "_", 
                       names = c("Condition", "Donor")) |> 
  # 1-3 are normal and 4-7 are oa
  mutate(Condition = ifelse(Condition %in% 1:3, "normal", "oa")) |> 
  mutate(Condition = factor(Condition, levels = c("normal", "oa")))

GSE168505_countMatrix <- GSE168505_raw_counts |> 
  column_to_rownames(var = "gene") |> 
  as.matrix() 
GSE168505_countMatrix_rownames <- rownames(GSE168505_countMatrix)
# Convert this matrix to integers
GSE168505_countMatrix <- apply(GSE168505_countMatrix, 2, as.integer)
rownames(GSE168505_countMatrix) <- GSE168505_countMatrix_rownames

# Make DESeqDataSet
GSE168505_dds <- DESeqDataSetFromMatrix(GSE168505_countMatrix, 
                                        GSE168505_colData, design = ~Condition) 

# Filter out lowly expressed genes
keep <- rowSums(counts(GSE168505_dds) >= 10) >= 1
GSE168505_dds <- GSE168505_dds[keep,]

# Fit model
GSE168505_dds <- DESeq(GSE168505_dds)
save(GSE168505_dds, file = "data/GSE168505/GSE168505_dds.rda")

# Shrink l2fc
GSE168505_shrink <- lfcShrink(GSE168505_dds,
                              coef = "Condition_oa_vs_normal") |>
  as_tibble(rownames = "symbol")

# All gene results
write_csv(GSE168505_shrink, file = "data/GSE168505/GSE168505_deseq_res.csv")

# Get significant genes; this study used padj < 0.05 and abs(log2FoldChange) > 2 significance cutoffs
GSE168505_pval05_l2fc2 <- GSE168505_shrink |> 
  filter(padj < 0.05 & abs(log2FoldChange) > 2) |> 
  write_csv("data/GSE168505/GSE168505_pval05_l2fc2.csv")

# Gene intersect
gene_intersect <- intersect(rownames(GSE114007_countMatrix),
                                    rownames(GSE168505_countMatrix))

GSE114007_GSE168404_rawCounts <- bind_cols(GSE114007_countMatrix[gene_intersect,], 
                                           GSE168505_countMatrix[gene_intersect,])
pcs_GSE114007_GSE168404 <- prcomp(t(GSE114007_GSE168404_rawCounts))

pc_df <- data.frame("PC1" = pcs_GSE114007_GSE168404$x[,1],
                    "PC2" = pcs_GSE114007_GSE168404$x[,2]) |> 
  rownames_to_column(var = "Sample") |> 
  mutate(Study = ifelse(Sample %in% c("1_S7",
                                      "2_S8",
                                      "3_S9",
                                      "4_S10",
                                      "5_S11",
                                      "6_S12",
                                      "7_S13"), "Fu", "Fisch"))


ggplot(data = pc_df, aes(x = PC1, 
                         y = PC2, col = Study)) +
  geom_point() 
