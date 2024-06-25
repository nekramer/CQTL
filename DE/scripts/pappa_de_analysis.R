library(tximeta)
# dbpyr version 2.3.4
library(DESeq2)
library(tidyverse)
library(ggrepel)

pappa_samplesheet <- read_csv("data/pappa_de/PAPPAsamplesheet.csv") |> 
  # Add Names and quant files
  mutate(names = paste0(Genotype, "_", Treatment, "_", Time),
         files = paste0("/work/users/n/e/nekramer/PAPPA_rna/bagPipes/output/quant/CQTL_CHON_",
                        names, "/quant.sf")) |> 
  dplyr::rename(Donor = Genotype) |> 
  mutate(Treatment = factor(Treatment, levels = c("vehicle", "PAPPA")))
# Import salmon transcript quantification-------------------------------------
se <- tximeta(pappa_samplesheet)

# Convert to gene-level scaled transcripts -------------------------------------
pappa_gse <- summarizeToGene(se)

save(pappa_gse, file = "data/pappa_de/pappa_gse.rda")
# Filter out lowly expressed genes
keep <- rowSums(assay(pappa_gse) >= 10) >= 2
pappa_gse_filtered <- pappa_gse[keep,]



# 3hr ---------------------------------------------------------------------
pappa_gse_3hr <- pappa_gse_filtered[, pappa_gse_filtered$Time == "3h"]
pappa_dds_3hr <- DESeqDataSet(pappa_gse_3hr, design = ~Donor + Treatment)

# Fit model
pappa_dds_3hr <- DESeq(pappa_dds_3hr)

# l2fc shrink
pappa_3hr_lfcshrink <- lfcShrink(pappa_dds_3hr, coef = "Treatment_PAPPA_vs_vehicle", format = "GRanges") |>
  plyranges::names_to_column("gene_id")

pappa_3hr_lfcshrink_normal <- lfcShrink(pappa_dds_3hr, 
                                        coef = "Treatment_PAPPA_vs_vehicle",
                                        type = "normal", format = "GRanges") |>
  plyranges::names_to_column("gene_id")

write_csv(as.data.frame(pappa_3hr_lfcshrink), file = "data/pappa_de/pappa_3hr_lfcshrink.csv")

pappa_3hr_p01_l2fc1 <- pappa_3hr_lfcshrink |> 
  filter(padj < 0.01 & abs(log2FoldChange) > 1)

pappa_3hr_lfcshrink <- pappa_3hr_lfcshrink |>
  as.data.frame() |>
  dplyr::rename(`3hr_apeglm_log2FoldChange` = log2FoldChange,
         `3hr_apeglm_lfcSE` = lfcSE,
         `3hr_pvalue` = pvalue,
         `3hr_padj` = padj) |> 
  na.omit()

pappa_3hr_lfcshrink_normal <- pappa_3hr_lfcshrink_normal |>
  as.data.frame() |>
  dplyr::select(gene_id, log2FoldChange, lfcSE) |> 
  dplyr::rename(`3hr_normal_log2FoldChange` = log2FoldChange,
         `3hr_normal_lfcSE` = lfcSE) |> 
  na.omit()

pappa_3hr <- full_join(pappa_3hr_lfcshrink, pappa_3hr_lfcshrink_normal, by = "gene_id") |> 
  dplyr::select(-baseMean)
# 6hr --------------------------------------------------------------------
pappa_gse_6hr <- pappa_gse_filtered[, pappa_gse_filtered$Time == "6h"]
pappa_dds_6hr <- DESeqDataSet(pappa_gse_6hr, design = ~Donor + Treatment)

# Fit model
pappa_dds_6hr <- DESeq(pappa_dds_6hr)

# l2fc shrink
pappa_6hr_lfcshrink <- lfcShrink(pappa_dds_6hr, coef = "Treatment_PAPPA_vs_vehicle", format = "GRanges") |>
  plyranges::names_to_column("gene_id")
pappa_6hr_lfcshrink_normal <- lfcShrink(pappa_dds_6hr, 
                                        coef = "Treatment_PAPPA_vs_vehicle",
                                        type = "normal", format = "GRanges") |>
  plyranges::names_to_column("gene_id")

write_csv(as.data.frame(pappa_6hr_lfcshrink), file = "data/pappa_de/pappa_6hr_lfcshrink.csv")

pappa_6hr_p01_l2fc1 <- pappa_6hr_lfcshrink |> 
  filter(padj < 0.01 & abs(log2FoldChange) > 1)

pappa_6hr_lfcshrink <- pappa_6hr_lfcshrink |>
  as.data.frame() |>
  dplyr::rename(`6hr_apeglm_log2FoldChange` = log2FoldChange,
                `6hr_apeglm_lfcSE` = lfcSE,
                `6hr_pvalue` = pvalue,
                `6hr_padj` = padj) |> 
  na.omit()

pappa_6hr_lfcshrink_normal <- pappa_6hr_lfcshrink_normal |>
  as.data.frame() |>
  dplyr::select(gene_id, log2FoldChange, lfcSE) |> 
  dplyr::rename(`6hr_normal_log2FoldChange` = log2FoldChange,
                `6hr_normal_lfcSE` = lfcSE) |> 
  na.omit()

pappa_6hr <- full_join(pappa_6hr_lfcshrink, pappa_6hr_lfcshrink_normal,
                       by = "gene_id") |> 
  dplyr::select(-baseMean)
# 24hr --------------------------------------------------------------------

pappa_gse_24hr <- pappa_gse_filtered[, pappa_gse_filtered$Time == "24h"]
pappa_dds_24hr <- DESeqDataSet(pappa_gse_24hr, design = ~Donor + Treatment)

# Fit model
pappa_dds_24hr <- DESeq(pappa_dds_24hr)

# l2fc shrink
pappa_24hr_lfcshrink <- lfcShrink(pappa_dds_24hr, coef = "Treatment_PAPPA_vs_vehicle", format = "GRanges") |>
  plyranges::names_to_column("gene_id")

pappa_24hr_lfcshrink_normal <- lfcShrink(pappa_dds_24hr, 
                                        coef = "Treatment_PAPPA_vs_vehicle",
                                        type = "normal", format = "GRanges") |>
  plyranges::names_to_column("gene_id")
write_csv(as.data.frame(pappa_24hr_lfcshrink), file = "data/pappa_de/pappa_24hr_lfcshrink.csv")

pappa_24hr_p01_l2fc1 <- pappa_24hr_lfcshrink |> 
  filter(padj < 0.01 & abs(log2FoldChange) > 1)

pappa_24hr_lfcshrink <- pappa_24hr_lfcshrink |>
  as.data.frame() |>
  dplyr::rename(`24hr_apeglm_log2FoldChange` = log2FoldChange,
                `24hr_apeglm_lfcSE` = lfcSE,
                `24hr_pvalue` = pvalue,
                `24hr_padj` = padj) |> 
  na.omit()

pappa_24hr_lfcshrink_normal <- pappa_24hr_lfcshrink_normal |>
  as.data.frame() |>
  dplyr::select(gene_id, log2FoldChange, lfcSE) |> 
  dplyr::rename(`24hr_normal_log2FoldChange` = log2FoldChange,
                `24hr_normal_lfcSE` = lfcSE) |> 
  na.omit()

pappa_24hr <- full_join(pappa_24hr_lfcshrink, pappa_24hr_lfcshrink_normal) |> 
  dplyr::select(-baseMean)




# Join all results --------------------------------------------------------

pappa_all <- full_join(pappa_3hr, pappa_6hr) |> 
  full_join(pappa_24hr) |> 
  dplyr::select(gene_id, seqnames, start, end, strand,
                `3hr_pvalue`, `3hr_padj`, `6hr_pvalue`, `6hr_padj`,
                `24hr_pvalue`, `24hr_padj`,
                `3hr_apeglm_log2FoldChange`, `3hr_apeglm_lfcSE`,
                `6hr_apeglm_log2FoldChange`, `6hr_apeglm_lfcSE`,
                `24hr_apeglm_log2FoldChange`, `24hr_apeglm_lfcSE`,
                `3hr_normal_log2FoldChange`, `3hr_normal_lfcSE`,
                `6hr_normal_log2FoldChange`, `6hr_normal_lfcSE`,
                `24hr_normal_log2FoldChange`, `24hr_normal_lfcSE`) |> 
  # Join back gene names from gse 
left_join(as.data.frame(rowData(pappa_gse_filtered)) |> dplyr::select("gene_id", "gene_name"),
          by = "gene_id") |> 
  relocate(gene_name, .after = "gene_id") |> 
  dplyr::rename(chrom = seqnames) |> 
  arrange(`3hr_padj`, `6hr_padj`, `24hr_padj`)

write_csv(pappa_all, file = "data/pappa_de/pappa_de_geneResults.csv")

# MA plot for all timepoints ----------------------------------------------

pappa_all_ma <- bind_rows(na.omit(pappa_3hr_ma), na.omit(pappa_6hr_ma), na.omit(pappa_24hr_ma)) |> 
  mutate(Time = factor(Time, levels = c("3hr", "6hr", "24hr"))) |> 
  # Join back gene names from gse 
  left_join(as.data.frame(rowData(pappa_gse_filtered)) |> dplyr::select("gene_id", "gene_name"),
            by = "gene_id") |> 
  mutate(dir_effect = ifelse(log2FoldChange > 0, "up", ifelse(log2FoldChange < 0, "down", "none"))) |> 
  mutate(sig_group = case_when(padj < 0.05 & dir_effect == "down" ~ "down_sig",
                               padj < 0.05 & dir_effect == "up" ~ "up_sig")) |> 
  mutate(sig_group = factor(sig_group))

ggplot(pappa_all_ma, mapping = aes(x = baseMean, y = log2FoldChange, col = sig_group)) +
  geom_point() +
  facet_wrap(vars(Time)) +
  scale_x_continuous(name = "mean of normalized counts",
                      trans = "log10") +
  scale_y_continuous(name = "shrunken (type = normal) log2FoldChange",
                     breaks = seq(-0.4, 0.4, 0.2)) +
  scale_color_manual(values = c("#78A1Cd", "#FBBE67")) +
  theme_minimal() +
  theme(text = element_text(family = "Helvetica"),
        legend.position = "none")

ggplot(pappa_all_ma |> filter(dir_effect != "none"), 
       mapping = aes(x = dir_effect, fill = dir_effect)) +
  geom_bar() +
  facet_wrap(vars(Time)) +
  scale_y_continuous(breaks = seq(0, 60, 10)) + 
  scale_fill_manual(values = c("#78A1Cd", "#FBBE67")) +
  theme_minimal()  +
  theme(text = element_text(family = "Helvetica"),
        legend.position = "none",
        axis.title.x = element_blank())
