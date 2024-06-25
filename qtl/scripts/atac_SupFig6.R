library(tidyverse)
library(data.table)
library(plyranges)
library(ggstar)
source("../plotting_utils.R")

# Overlapping peaks with eQTLs --------------------------------------------

# Get significant differential ATAC peaks (abs log2FC > 1 and padj < 0.01)
load("data/atac/atac_res.rda")
diff_atac_peaks <- atac_res |> 
  filter(abs(log2FoldChange) > 1 & padj < 0.01)

# Get candidate variants for each independent eQTL signal, defined as
# variants that pass their associated eGene's nominal p-value threshold
PBS_FNF_signals_eSNPs <- list()

for (chrom in 1:22){
  
  # PBS
  pbs_chrom_signal_eSNPs <- fread(paste0("data/eqtl/qtl_nom/CTL_PEER_k20_genoPC_allSignals_nom1Mb_MAFs_chr",
                                         chrom, ".csv"), data.table = FALSE) |> 
    filter(nom_sig == 1) |> 
    mutate(condition = "PBS")
  
  # FN-f
  fnf_chrom_signal_eSNPs <- fread(paste0("data/eqtl/qtl_nom/FNF_PEER_k22_genoPC_allSignals_nom1Mb_MAFs_chr",
                                         chrom, ".csv"), data.table = FALSE) |> 
    filter(nom_sig == 1) |> 
    mutate(condition = "FN-f")
  
  # Combine and append to list
  PBS_FNF_signals_eSNPs[[paste0("chr", chrom)]] <- bind_rows(pbs_chrom_signal_eSNPs,
                                                             fnf_chrom_signal_eSNPs)
}

PBS_FNF_signals_eSNPs <- bind_rows(PBS_FNF_signals_eSNPs) |> 
  distinct(gene_id, variantID, .keep_all = TRUE)

# Create GRanges for overlapping
PBS_FNF_signals_eSNPs_gr <- makeGRangesFromDataFrame(PBS_FNF_signals_eSNPs, keep.extra.columns = TRUE,
                                                     seqnames.field = "variant_chr", start.field = "variant_start",
                                                     end.field = "variant_end")


# Find overlaps between eSNPs and all peaks
eSNP_peak_overlaps <- findOverlaps(PBS_FNF_signals_eSNPs_gr, atac_res)
eSNP_peaks <- atac_res[unique(subjectHits(eSNP_peak_overlaps))] |> 
  as_tibble()


# Number and percentage of unique eSNPs overlapping any peak
no_eSNP_peak_overlaps <- queryHits(eSNP_peak_overlaps) |> unique () |> length()
perc_eSNP_peak_overlaps <- (no_eSNP_peak_overlaps/length(PBS_FNF_signals_eSNPs_gr))*100

# Find overlaps between eSNPs and differential peaks
eSNP_diff_peak_overlaps <- findOverlaps(PBS_FNF_signals_eSNPs_gr, diff_atac_peaks)

# Number and percentage of unique eSNPs overlapping a differential peak
no_eSNP_diff_peak_overlaps <- queryHits(eSNP_diff_peak_overlaps) |> unique() |> length()
perc_eSNP_diff_peak_overlaps <- no_eSNP_diff_peak_overlaps/length(PBS_FNF_signals_eSNPs_gr)*100


# Peak overlap with PBS-specific, FN-f-specific, and shared lead eQTLs and  
# LD > 0.8 --------

# Get all ATAC peaks
load("data/atac/atac_res.rda")

#### Creating background: any lead eQTL and LD > 0.8 that overlaps an ATAC peak

# PBS eQTLs and LD > 0.8
pbs_eqtls_LD <- fread("data/eqtl/CTL_PEER_k20_genoPC_cond1Mb_topSignals_rsID_LD.csv",
                      data.table = FALSE) |> 
  filter(R2 > 0.8) |> 
  separate_wider_delim(cols = "ld_variantID", delim = ":", names = c(NA, "ld_pos", NA, NA), cols_remove = FALSE)

# FNF eQTLs and LD > 0.8
fnf_eqtls_LD <- fread("data/eqtl/FNF_PEER_k22_genoPC_cond1Mb_topSignals_rsID_LD.csv",
                      data.table = FALSE) |> 
  filter(R2 > 0.8) |> 
  separate_wider_delim(cols = "ld_variantID", delim = ":", names = c(NA, "ld_pos", NA, NA), cols_remove = FALSE)

# Bind PBS and FN-f
PBS_FNF_eqtls_LD <- bind_rows(pbs_eqtls_LD, fnf_eqtls_LD) |> 
  distinct(gene_id, ld_variantID, .keep_all = TRUE)

# Create GRanges for overlapping
PBS_FNF_eqtls_LD_gr <- makeGRangesFromDataFrame(PBS_FNF_eqtls_LD, keep.extra.columns = TRUE,
                                                     seqnames.field = "variant_chr", start.field = "ld_pos",
                                                     end.field = "ld_pos")


# Find overlaps between eSNPs and all peaks
eqtls_LD_peak_overlaps <- findOverlaps(PBS_FNF_eqtls_LD_gr, atac_res)
eqtl_LD_peaks <- atac_res[unique(subjectHits(eqtls_LD_peak_overlaps))] |> 
  as_tibble()

#### Target: Peaks overlapping high confidence PBS-specific, FN-f-specific,
# and shared eQTLs and LD > 0.8

# Get high confidence eGenes
pbs_highconf_regenes <- read_csv("data/reqtl/CTL_sig01_beta_donor_reQTLs_PEER_k20_genoPC.csv")
fnf_highconf_regenes <- read_csv("data/reqtl/FNF_sig01_beta_donor_reQTLs_PEER_k22_genoPC.csv")

# Shared eGenes
pbs_egenes <- read_csv("data/eqtl/CTL_PEER_k20_genoPC_perm1Mb_sig_rsID.csv") |> 
  pull(gene_id)
fnf_egenes <- read_csv("data/eqtl/FNF_PEER_k22_genoPC_perm1Mb_sig_rsID.csv") |> 
  pull(gene_id)

shared_egenes <- intersect(pbs_egenes, fnf_egenes)

# High confidence PBS and LD > 0.8
highconf_PBS_signalSNPs_LD <- fread("data/eqtl/CTL_PEER_k20_genoPC_cond1Mb_topSignals_rsID_LD.csv", 
                                 data.table = FALSE) |> 
  filter(gene_id %in% pbs_highconf_regenes$gene_id) |> 
  filter(R2 > 0.8) |> 
  separate_wider_delim(cols = "ld_variantID", delim = ":", names = c(NA, "ld_pos", NA, NA), cols_remove = FALSE)

# High confidence FNF and LD > 0.8
highconf_FNF_signalSNPs_LD <- fread("data/eqtl/FNF_PEER_k22_genoPC_cond1Mb_topSignals_rsID_LD.csv",
                                    data.table = FALSE) |> 
  filter(gene_id %in% fnf_highconf_regenes$gene_id) |> 
  filter(R2 > 0.8) |> 
  separate_wider_delim(cols = "ld_variantID", delim = ":", names = c(NA, "ld_pos", NA, NA), cols_remove = FALSE)


# Shared from either PBS or FNF
shared_signalSNPs_LD <- bind_rows(fread("data/eqtl/CTL_PEER_k20_genoPC_cond1Mb_topSignals_rsID_LD.csv", 
                                             data.table = FALSE) |> 
                                         filter(gene_id %in% shared_egenes) |> 
                                         filter(R2 > 0.8) |> 
                                         separate_wider_delim(cols = "ld_variantID", 
                                                              delim = ":",
                                                              names = c(NA, "ld_pos", NA, NA), cols_remove = FALSE) |> 
                                         dplyr::select(gene_id, ld_variantID, gene_chr, ld_pos),
                                       fread("data/eqtl/FNF_PEER_k22_genoPC_cond1Mb_topSignals_rsID_LD.csv", 
                                             data.table = FALSE) |> 
                                         filter(gene_id %in% shared_egenes) |>
                                         filter(R2 > 0.8) |> 
                                         separate_wider_delim(cols = "ld_variantID", 
                                                              delim = ":",
                                                              names = c(NA, "ld_pos", NA, NA), cols_remove = FALSE) |>
                                         dplyr::select(gene_id, ld_variantID, gene_chr, ld_pos)) |> 
  distinct()

# Make GRanges for overlapping
highconf_PBS_signalSNPs_LD_gr <- makeGRangesFromDataFrame(highconf_PBS_signalSNPs_LD,
                                                          keep.extra.columns = TRUE,
                                                          seqnames.field = "variant_chr", 
                                                          start.field = "ld_pos",
                                                          end.field = "ld_pos")
highconf_FNF_signalSNPs_LD_gr <- makeGRangesFromDataFrame(highconf_FNF_signalSNPs_LD,
                                                          keep.extra.columns = TRUE,
                                                          seqnames.field = "variant_chr", 
                                                          start.field = "ld_pos",
                                                          end.field = "ld_pos")
shared_signalSNPs_LD_gr <- makeGRangesFromDataFrame(shared_signalSNPs_LD,
                                                          keep.extra.columns = TRUE,
                                                          seqnames.field = "gene_chr", 
                                                          start.field = "ld_pos",
                                                          end.field = "ld_pos")

# Find overlaps between these eSNPs and all peaks
highconf_PBS_LD_peak_overlaps <- findOverlaps(highconf_PBS_signalSNPs_LD_gr, atac_res)
highconf_FNF_LD_peak_overlaps <- findOverlaps(highconf_FNF_signalSNPs_LD_gr, atac_res)
shared_LD_peak_overlaps <- findOverlaps(shared_signalSNPs_LD_gr, atac_res)

# Subset atac results for overlaps
highconf_PBS_LD_peaks <- atac_res[unique(subjectHits(highconf_PBS_LD_peak_overlaps))] |> 
  mutate(condition = "PBS") |> 
  mutate(diff_peak = ifelse(padj < 0.01 & abs(log2FoldChange) > 1, "diff", "not_diff")) |> 
  as_tibble()

highconf_FNF_LD_peaks <- atac_res[unique(subjectHits(highconf_FNF_LD_peak_overlaps))] |> 
  mutate(condition = "FN-f") |> 
  mutate(diff_peak = ifelse(padj < 0.01 & abs(log2FoldChange) > 1, "diff", "not_diff")) |> 
  as_tibble()

shared_LD_peaks <- atac_res[unique(subjectHits(shared_LD_peak_overlaps))] |> 
  mutate(condition = "shared") |> 
  mutate(diff_peak = ifelse(padj < 0.01 & abs(log2FoldChange) > 1, "diff", "not_diff")) |> 
  as_tibble()

# Wilcox testing 
# PBS-specific enriched for less accessible peaks compared to all peaks
# # that overlap any eQTL?
eqtl_LD_notPBS <- anti_join(eqtl_LD_peaks,
                               highconf_PBS_LD_peaks |>
                                 dplyr::select(-condition,-diff_peak))

pbs_ld_peak_wilcox <- wilcox.test(x = highconf_PBS_LD_peaks$log2FoldChange,
                               y = eqtl_LD_notPBS$log2FoldChange,
                               alternative = "less")

# # FN-f-specific enriched for more accessible peaks compared to all peaks
# # that overlap any eQTL?
eqtl_LD_notFNF <- anti_join(eqtl_LD_peaks,
                               highconf_FNF_LD_peaks |>
                                 dplyr::select(-condition,-diff_peak))

fnf_ld_peak_wilcox <- wilcox.test(x = highconf_FNF_LD_peaks$log2FoldChange,
                               y = eqtl_LD_notFNF$log2FoldChange,
                               alternative = "greater")

# Join all condition groups for plotting
PBS_FNF_shared_LD_peaks <- bind_rows(highconf_PBS_LD_peaks,
                                     highconf_FNF_LD_peaks,
                                     shared_LD_peaks) |> 
  mutate(condition = factor(condition, levels = c("PBS", "shared", "FN-f")))

highconf_LD_peakOverlap_boxplots <- ggplot(PBS_FNF_shared_LD_peaks,
       aes(x = condition, y = log2FoldChange)) +
  geom_hline(yintercept = 0, lty = 2, linewidth = 0.25) +
  geom_jitter(aes(color = condition), position = position_jitter(width = 0.1), size = 0.25,
              alpha = 0.3, pch = 16) +
  geom_boxplot(aes(color = condition, fill = condition), outlier.shape = NA, alpha = 0.8) +

  scale_color_manual(values = c("#48617b", "#746541", "#97723e")) +
  scale_x_discrete(labels = c("CA regions<br> that overlap<br> **PBS-specific**<br> eQTLs",
                              "CA regions<br> that overlap<br> **shared**<br> eQTLs",
                              "CA regions<br> that overlap<br> **FN-f-specific**<br> eQTLs")) +
  scale_fill_manual(values = c("#78A1Cd", "#E1D8C7", "#FBBE67"))  + 
  scale_y_continuous(name = "log~2~ Fold Change", expand = c(0, 0),
                     limits = c(-3, 5.5)) +
  annotate("text",
           label = paste0("pval = ", format(pbs_ld_peak_wilcox$p.value, digits = 3)),
           x = 1, y = Inf, family = "Helvetica", size = 2.25, vjust = 1) +
  annotate("text",
           label = paste0("pval = ",format(fnf_ld_peak_wilcox$p.value, digits = 3)),
           x = 3, y = Inf, family = "Helvetica", size = 2.25, vjust = 1) +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        axis.ticks.x = element_blank(),
        axis.line.y = element_line(color = "black", 
                                   linewidth = 0.2),
        axis.ticks.length.y = unit(-0.1, "cm"),
        axis.ticks.y = element_line(color = "black",
                                    linewidth = 0.2),
        panel.background = element_rect(fill = 'transparent', color = "transparent"),
        plot.background = element_rect(fill = 'transparent', color = "transparent"),
        axis.text = element_text(color = "black"),
        axis.text.x = element_markdown(margin = margin(t = 3), family = "Helvetica", size = 7),
        text = element_text(family = "Helvetica", size = 8),
        axis.title.y = element_markdown(size = 7, margin = margin(r = -2))) +
  coord_cartesian(clip = "off")

save(highconf_LD_peakOverlap_boxplots, 
     file = "plots/atac/highconf_LD_peakOverlap_boxplots.rda")
# Create BED files for TF motif enrichment with Homer ---------------------

# Target: peaks that overlap PBS-specific lead eQTLs and LD > 0.8
highconf_PBS_LD_peaks |> 
  dplyr::select(seqnames, start, end, strand) |> 
  mutate(peakID = paste0(seqnames, "_", start, "_", end),
         dummy = NA) |> 
  relocate(peakID, .after = end) |> 
  relocate(dummy, .before = strand) |> 
  write_delim(file = "data/atac/homer/peaks_highconf_LD_PBS.txt",
              delim = "\t", col_names = FALSE)


# Target: peaks that overlap FNF-specific lead eQTLs and LD > 0.8
highconf_FNF_LD_peaks |> 
  dplyr::select(seqnames, start, end, strand) |> 
  mutate(peakID = paste0(seqnames, "_", start, "_", end),
         dummy = NA) |> 
  relocate(peakID, .after = end) |> 
  relocate(dummy, .before = strand) |> 
  write_delim(file = "data/atac/homer/peaks_highconf_LD_FNF.txt",
              delim = "\t", col_names = FALSE)

# Background : all peaks that overlap any lead eQTL w/ LD > 0.8
eqtl_LD_peaks |> 
  as_tibble() |> 
  dplyr::select(seqnames, start, end, strand) |> 
  mutate(peakID = paste0(seqnames, "_", start, "_", end),
         dummy = NA) |> 
  relocate(peakID, .after = end) |> 
  relocate(dummy, .before = strand) |> 
  write_delim(file = "data/atac/homer/peaks_all_eQTL_LD.txt",
              delim = "\t", col_names = FALSE)

# Launch homer, PBS LD w/ peaks that overlap any eQTL and LD
system("scripts/peaks_eqtl_homer.sh data/atac/homer/peaks_highconf_LD_PBS.txt data/atac/homer/peaks_all_eQTL_LD.txt data/atac/homer/PBS_LD")

# Launch homer, FNF LD w/ peaks that overlap any eQTL and LD
system("scripts/peaks_eqtl_homer.sh data/atac/homer/peaks_highconf_LD_FNF.txt data/atac/homer/peaks_all_eQTL_LD.txt data/atac/homer/FNF_LD")

# Barplot of peak TF motifs -----------------------------------------------

# Peaks that overlap PBS-specific eQTLs
pbs_peak_motifs <- read_delim("data/atac/homer/PBS_LD/knownResults.txt") |> 
  # Convert percentages to numbers
  mutate(across(c(`% of Target Sequences with Motif`, 
                  `% of Background Sequences with Motif`),
                ~ gsub("%", "", .))) |> 
  mutate(across(c(`% of Target Sequences with Motif`, 
                  `% of Background Sequences with Motif`), as.numeric)) |> 
  # Calculate log2 enrichment
  mutate(log2enrichment = log2(`% of Target Sequences with Motif`/`% of Background Sequences with Motif`)) |> 
  # Calculate -log10pval
  mutate(log10pval = -log10(exp(`Log P-value`))) |> 
  slice_max(order_by = log10pval, n = 6) |> 
  mutate(condition = "PBS")

# Remove second p53
pbs_peak_motifs <- pbs_peak_motifs[-5,]
# Pull out first part of motif name and convert to all uppercase
pbs_peak_motifs$Name <- toupper(unlist(lapply(str_split(pbs_peak_motifs$`Motif Name`, 
                                                  '[(]'), `[[`, 1)))



# Peaks that overlap FNF-specific eQTLs
fnf_peak_motifs <- read_delim("data/atac/homer/FNF_LD/knownResults.txt") |> 
  # Convert percentages to numbers
  mutate(across(c(`% of Target Sequences with Motif`, 
                  `% of Background Sequences with Motif`),
                ~ gsub("%", "", .))) |> 
  mutate(across(c(`% of Target Sequences with Motif`, 
                  `% of Background Sequences with Motif`), as.numeric)) |> 
  # Calculate log2 enrichment
  mutate(log2enrichment = log2(`% of Target Sequences with Motif`/`% of Background Sequences with Motif`)) |> 
  # Calculate -log10pval
  mutate(log10pval = -log10(exp(`Log P-value`))) |> 
  slice_max(order_by = log10pval, n = 5) |> 
  mutate(condition = "FN-f")

# Pull out first part of motif name and convert to all uppercase
fnf_peak_motifs$Name <- toupper(unlist(lapply(str_split(fnf_peak_motifs$`Motif Name`, 
                                                        '[(]'), `[[`, 1)))

# Join PBS and FNF for plotting
pbs_fnf_peak_motifs <- bind_rows(pbs_peak_motifs, fnf_peak_motifs) |> 
  arrange(log10pval) |> 
  mutate(Name = factor(Name, levels = Name)) |> 
  mutate(condition = factor(condition, levels = c("PBS", "FN-f")))

# Barplot
highconf_peakOverlap_tfmotifs <- ggplot(pbs_fnf_peak_motifs, aes(x = log10pval, y = Name, fill = condition)) +
  geom_bar(stat = "identity") +
  facet_wrap(~condition, ncol = 1, scales = "free_y",
             labeller = as_labeller(c(`PBS` = "Chromatin accessible peaks overlapping\n PBS-specific eQTLs",
                                      `FN-f` = "Chromatin accessible peaks overlapping\n FN-f-specific eQTLs"))) +
  geom_text(aes(x = 0, label = Name), hjust = 0, family = "Helvetica",
            size = 3) +
  scale_fill_manual(values = c("#78A1Cd", "#FBBE67")) +
  scale_x_continuous(name = "-log~10~pval", expand = c(0, 0),
                     limits = c(0, 4)) +
  theme(panel.background = element_rect(fill = 'transparent', color = "transparent"),
        plot.background = element_rect(fill = 'transparent', color = "transparent"),
        text = element_text(family = "Helvetica"),
        legend.position = "None",
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_markdown(size = 6, margin = margin(t = 0)),
        axis.text.x = element_text(color = "black", size = 5, margin = margin(t = 0)),
        strip.background = element_blank(),
        strip.placement = "outside",
        axis.ticks = element_blank(),
        axis.line.x = element_line(linewidth = 0.25),
        strip.text = element_text(size = 6, color = "black", face = "bold", hjust = 0.5,
                                  margin = margin(b = -1)),
        panel.spacing = unit(0, "mm"), 
        plot.title = element_text(hjust = 0.5, face = "bold", size = 7, margin = margin(b = 0.5))) +
  coord_cartesian(clip = "off") +
  ggtitle("Motif Enrichment")

save(highconf_peakOverlap_tfmotifs, file = "plots/atac/highconf_peakOverlap_tfmotifs.rda")


# plotgardener layout -----------------------------------------------------

pdf(file = "plots/atac/SupFig5.pdf", width = 5.25, height = 3.5)
pageCreate(width = 5.25, height = 3.5, showGuides = FALSE)
plotText("A", x = 0.1, y = 0.1, just = c("left", "top"), fontfamily = "Helvetica",
         fontsize = 14, fontface = "bold")
plotGG(highconf_LD_peakOverlap_boxplots, x = 0.1, y = 0.1, width = 2.75, height = 3.2)

plotText("B", x = 2.9, y = 0.1, just = c("left", "top"), fontfamily = "Helvetica",
         fontsize = 14, fontface = "bold")
plotGG(highconf_peakOverlap_tfmotifs, x = 3, y = 0, width = 2.2, height = 3.4)
dev.off()


# Differential ATAC peaks table -------------------------------------------

# Get significant differential ATAC peaks (abs log2FC > 1 and padj < 0.01)
load("data/atac/atac_res.rda")
diff_atac_peaks <- atac_res |> 
  filter(abs(log2FoldChange) > 1 & padj < 0.01) |> 
  as_tibble() |> 
  dplyr::select(-strand, -pvalue) |> 
  dplyr::rename(chrom = seqnames)

write_csv(diff_atac_peaks, file = "tables/SupTable11.csv")

ss <- gs4_create(name = "SupTable11_differential_ATAC_peaks")
write_sheet(diff_atac_peaks,
            ss, sheet = "FNF_vs_PBS")

drive_mv(file = "SupTable11_differential_ATAC_peaks", 
         path = as_dribble("chQTL/chQTL paper/Figures and Tables"))

