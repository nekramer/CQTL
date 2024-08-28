library(tidyverse)
library(data.table)
library(plotgardener)
library(scales)
library(GenomicRanges)
library(mariner)
library(InteractionSet)
library(nullranges)
library(org.Hs.eg.db)
library(ggtext)
library(stats)
source("../plotting_utils.R")
source("../utils.R")

# Functions ---------------------------------------------------------------

# Function to find any overlapping loops with a gene based on promoters
# @param gene Gene symbol
# @param loops GInteractions of loop list to search
# @param loop_binsize The binsize of loop anchors to use when checking connections, i.e. expanding anchors to 30Kb
# @param promoters_data GRanges of gene promoters data to use for looking up gene and then connecting to loops
# @param txdb TxDb object to use for getting a gene's transcripts
get_gene_loops <- function(gene, loops, loop_binsize, promoters_data, 
                           txdb = ensembl_txdb){
  
  # Get all the transcripts for that gene (calling function to get ENSEMBL ID)
  gene_transcripts <- AnnotationDbi::select(txdb,
                                           keys = symbol_to_ensembl(gene), 
                                           keytype = "GENEID",
                                           columns = "TXID")
  
  # Filter promoters for gene promoters based on gene_transcripts
  gene_promoters <- promoters_data |> 
    filter(tx_id %in% gene_transcripts$TXID)
  
  # Adjust binsize of loops
  binsize_loops <- snapToBins(loops, binSize = loop_binsize)
  
  
  # Get overlaps between the gene promoters and loops
  gene_loop_overlaps <- findOverlaps(query = binsize_loops,
                                    subject = gene_promoters)
  
  # Subset loop data for these overlaps
  gene_loops <- loops[unique(queryHits(gene_loop_overlaps))]
  return(gene_loops)
}

# Function to get Hi-C loops that connect an eQTL signal and its eGene
get_signal_eGene_loops <- function(eSNP_eGene, promoters_mapped_genes, loops, 
                                   loop_binsize){
  
  eGene_id <- eSNP_eGene[["gene_id"]]
  signal_start <- as.numeric(eSNP_eGene[["min_signal_pos"]])
  signal_end <- as.numeric(eSNP_eGene[["max_signal_pos"]])
  strand <- eSNP_eGene[["gene_strand"]]
  chrom <- eSNP_eGene[["gene_chr"]]
  
  # Define signal range
  signal_range <- GRanges(seqnames = chrom,
                          ranges = IRanges(start = signal_start,
                                           end = signal_end),
                          strand = strand)
  
  # Grab eGene_id promoter ranges
  eGene_promoters <- promoters_mapped_genes |> 
    filter(GENEID == eGene_id)
  
  # Get loop interactions that link signal and any of the eGene's promoters
  # linkOverlaps returns a DataFrame with the columns `query`, `subject1`, and `subject2`
  loop_interactions <- linkOverlaps(snapToBins(loops, binSize = loop_binsize), 
                                    subject1 = signal_range, 
                                    subject2 = eGene_promoters)
  
  if (nrow(loop_interactions) > 0){
    return(loops[unique(loop_interactions$query)])
  } else {
    return(NA)
  }
  
}


# RNA log2FC AT STATIC VS GAINED LOOP ANCHORS -----------------------------
CQTL_10kb_gained_static_loops_expression <- read_csv("/proj/phanstiel_lab/Data/processed/CQTL/hic/loops/CQTL_10kb_gained_static_loops_expression.csv")


## Density plot 
anchor_expression_density <- ggplot(CQTL_10kb_gained_static_loops_expression, 
                                    aes(x = log2FoldChange, group = loopCategory,
                                        fill = loopCategory)) +
  geom_vline(xintercept = 0, lty = 2, lwd = 0.25) +
  geom_density(alpha = 0.5, lwd = 0.25) +
  geom_richtext(data = data.frame(log2FoldChange = -2.5, 
                              y = 0.15,
                              label = 'genes at <span style="color:grey30">**static**</span><br>loop anchors'),
            aes(x = log2FoldChange, y = y, label = label), inherit.aes = FALSE,
            label.color = NA, family = "Helvetica", hjust = 1, fill = "transparent", 
            size = 2.75) +
  geom_richtext(data = data.frame(log2FoldChange = 3.5, 
                                  y = 0.15,
                                  label = 'genes at <span style="color:#267C80">**gained**</span><br>loop anchors'),
                aes(x = log2FoldChange, y = y, label = label), inherit.aes = FALSE,
                label.color = NA, family = "Helvetica", hjust = 0, fill = "transparent",
                size = 2.75) +
  scale_fill_manual(values = c("grey30", "#267C80")) +
  scale_y_continuous(expand = c(0, 0), name = "density") +
  scale_x_continuous(expand = c(0,0), 
                     limits = c(-11, 15), breaks = seq(-10, 15, 5)) +
  xlab("RNA log~2~FC") +
  theme_custom_scatterplot() +
  theme(legend.position = "none",
        axis.ticks.x = element_blank(),
        axis.text = element_text(color = "black", size = 6),
        axis.title.x = element_markdown(size = 7, margin = margin(t = -0.5)),
        axis.title.y = element_text(size = 7),
        plot.margin = margin(r = 15, t = 10, b = 10, l = 10))


save(anchor_expression_density, 
     file = "plots/distal_Fig4/anchor_expression_density.rda")

# GAINED/STATIC LOOP ANCHOR GENE GO TERMS -------------------------------------

# Get reduced, significant GO terms
gained_loop_go_data <- 
  read_delim("/proj/phanstiel_lab/Data/processed/CQTL/hic/homer/biological_process.txt") %>%
  mutate(pval = exp(1)^logP) %>%
  filter(pval < 0.01)
gained_loop_go <- reduceGO(gained_loop_go_data,
                     category = "gained")
  
# Grab terms for plotting  
gained_loop_go_plotting <- gained_loop_go |> 
  filter(parentTerm %in% c("regulation of inflammatory response",
                           "reactive oxygen species metabolic process",
                           "regulation of catabolic process",
                           "regulation of response to stress",
                           "extracellular matrix disassembly", 
                           "immune system process",
                           "interleukin-8 production",
                           "macrophage activation",
                           "regulation of CD40 signaling pathway",
                           "negative regulation of bone resorption")) |> 
  arrange(`-log10pval`)

# Set factors for ordering bars
gained_loop_go_plotting$parentTerm <- factor(gained_loop_go_plotting$parentTerm,
                                    levels = unique(gained_loop_go_plotting$parentTerm))

loop_go_barplot <- ggplot(gained_loop_go_plotting, aes(x = `-log10pval`, y = parentTerm)) +
  geom_bar(stat = "identity", fill = "#267C80", alpha = 0.5) +
  scale_x_continuous(limits = c(0, 4), expand = c(0, 0), name = "-log~10~pval",
                     breaks = seq(0, 4, 1)) +
  geom_text(aes(x = 0, label = parentTerm), hjust = 0, family = "Helvetica",
            size = 2.5) +
  theme(panel.background = element_rect(fill = 'transparent', color = "transparent"),
        plot.background = element_rect(fill = 'transparent', color = "transparent"),
        text = element_text(family = "Helvetica"),
        legend.position = "None",
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_markdown(size = 7, margin = margin(t = -3)),
        axis.text.x = element_text(color = "black", size = 6, margin = margin(t = -1)),
        axis.ticks = element_blank(),
        axis.line.x = element_line(linewidth = 0.25),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 8, margin = margin(b = -1))) +
  ggtitle("GO Terms")

save(loop_go_barplot, file = "plots/distal_Fig4/loop_go_barplot.rda")


# CONTACT FREQUENCY BETWEEN DISTAL LEAD ESNPS AND THEIR EGENE PROMOTER ---------

#distal_leadVar_pixel_matches <- load("data/hic/distal_leadVar_pixel_matches.rda")
load("data/hic/distal_leadVar_ld_pixel_matches.rda")

# Format data for boxplots   
distal_leadVar_ld_pixel_contactfreq <- c(focal(distal_leadVar_ld_pixel_matches), 
                                      distal_leadVar_ld_pixel_matches) |> 
  # Make tibble for tidy manipulation
  as_tibble() |> 
  # Grab relevant columns (ld_variantID, gene_id, gene_category, log2_contact_freq)
  dplyr::select(ld_variantID, gene_id, gene_category, mean_contact_freq) |>
  # Make group names based on gene_category
  mutate(gene_category_label = ifelse(gene_category == "eGene", "contact between<br> distal SNP<br> and eGene", 
                                    "contact between<br> distal SNP<br> and distance-matched<br> gene")) |> 
  # Change to factor correct plotting order
  mutate(gene_category_label = factor(gene_category_label, 
                                    levels = c("contact between<br> distal SNP<br> and eGene", 
                                               "contact between<br> distal SNP<br> and distance-matched<br> gene")))
# Calculate p-value
wilcox_contact_freq_test <- wilcox.test(x = distal_leadVar_ld_pixel_contactfreq |> 
                                          filter(gene_category == "eGene") |> 
                                          pull(mean_contact_freq),
                                        y = distal_leadVar_ld_pixel_contactfreq |>  
                                          filter(gene_category == "other") |> 
                                          pull(mean_contact_freq),
                                        alternative = "greater")

# Plot
distal_ld_pixel_contactfreq_boxplot <- ggplot(distal_leadVar_ld_pixel_contactfreq, 
                                             mapping = aes(x = gene_category_label, 
                                                           y = mean_contact_freq)) +
  geom_boxplot(aes(color = gene_category_label, fill = gene_category_label), outlier.shape = NA) +
  annotate("text",
           label = paste0("pval = ", format(wilcox_contact_freq_test$p.value, digits = 3)),
           x = 1.5, y = Inf, family = "Helvetica", size = 2, vjust = 1) +
  scale_fill_manual(values = c("#6DB2B6", "grey90")) +
  scale_color_manual(values = c("#32696C", "grey20")) +
  scale_y_continuous(name = "contact frequency", expand = c(0,0),
                     trans = scales::pseudo_log_trans(base = 2),
                     breaks = c(0, 5, 10, 20, 40, 80, 160, 320)) +
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
        axis.text.x = element_markdown(margin = margin(t = 3), family = "Helvetica", size = 6),
        text = element_text(family = "Helvetica", size = 8),
        axis.title.y = element_markdown(size = 7, margin = margin(r = -2))) +
  coord_cartesian(clip = "off")

save(distal_ld_pixel_contactfreq_boxplot, 
     file = "plots/distal_Fig4/distal_ld_pixel_contactfreq_boxplot.rda")  

# CONTACT FREQUENCY OF DISTAL PBS-SPECIFIC, SHARED, AND FNF-SPECIFIC EQTLS -----

eqtl_categories_distal_egene_pixels_counts <- read_csv("/proj/phanstiel_lab/Data/processed/CQTL/hic/eqtl_categories_distal_egene_pixels_counts.csv")

## Wilcox testing

# PBS-specific enriched for higher PBS contact freq compared to shared distal eqtls?
pbs_contactfreq_wilcox <- wilcox.test(x = eqtl_categories_distal_egene_pixels_counts |> 
                                        filter(eqtl_category == "PBS-specific") |> 
                                        pull(mean_log2FC_contact_freq),
                               y = eqtl_categories_distal_egene_pixels_counts |> 
                                 filter(eqtl_category == "shared") |> 
                                 pull(mean_log2FC_contact_freq),
                               alternative = "less")
# FNF-specific enriched for higher FNF contact freq compared to shared distal eqtls?
fnf_contactfreq_wilcox <- wilcox.test(x = eqtl_categories_distal_egene_pixels_counts |> 
                                        filter(eqtl_category == "FNF-specific") |> 
                                        pull(mean_log2FC_contact_freq),
                                      y = eqtl_categories_distal_egene_pixels_counts |> 
                                        filter(eqtl_category == "shared") |> 
                                        pull(mean_log2FC_contact_freq),
                                      alternative = "greater")

## Boxplot

eqtl_cat_distal_contactfreq_boxplot <- ggplot(eqtl_categories_distal_egene_pixels_counts,
       mapping = aes(x = eqtl_category_label,
                     y = mean_log2FC_contact_freq)) +
  geom_hline(yintercept = 0, lty = 2, linewidth = 0.25) +
  geom_boxplot(aes(color = eqtl_category_label, fill = eqtl_category_label), outlier.shape = NA,
               alpha = 1) +
  annotate("text",
           label = paste0("pval = ", format(pbs_contactfreq_wilcox$p.value, digits = 3)),
           x = 1, y = Inf, family = "Helvetica", size = 2, vjust = 1) +
  annotate("text",
           label = paste0("pval = ", format(fnf_contactfreq_wilcox$p.value, digits = 3)),
           x = 3, y = Inf, family = "Helvetica", size = 2, vjust = 1) +
  
  scale_fill_manual(values = c("#78A1Cd", "#E1D8C7", "#FBBE67")) +
  scale_color_manual(values = c("#48617b", "#746541", "#97723e")) +
  scale_y_continuous(name = "log~2~FC contact frequency", expand = c(0, 0), 
                     limits = c(-0.9, 0.9), breaks = c(-0.8, -0.6, -0.4, -0.2, 0,
                                                        0.2, 0.4, 0.6, 0.8)) +
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
        axis.text.x = element_markdown(margin = margin(t = 3), family = "Helvetica", size = 6),
        text = element_text(family = "Helvetica", size = 8),
        axis.title.y = element_markdown(size = 7, margin = margin(r = -2))) +
  coord_cartesian(clip = "off")

save(eqtl_cat_distal_contactfreq_boxplot, 
     file = "plots/distal_Fig4/eqtl_cat_distal_contactfreq_boxplot.rda")  
  
# GETTING LOOPS FOR JUN, IL6, AND MMP13 ----------------------------------

# Ensembl TxDb
ensembl_txdb <- makeTxDbFromEnsembl(organism = "Homo sapiens", release = 111)
seqlevelsStyle(ensembl_txdb) <- "UCSC"
load("/proj/phanstiel_lab/Data/processed/CQTL/hic/loops/CQTL_10kb_sig_gained_loops.rda")
ensembl_promoters <- promoters(ensembl_txdb)

jun_loops <- get_gene_loops(gene = "JUN", loops = CQTL_10kb_sig_gained_loops, 
                            loop_binsize = 30e3, promoters_data = ensembl_promoters)
il6_loops <- get_gene_loops(gene = "IL6", loops = CQTL_10kb_sig_gained_loops, 
                            loop_binsize = 30000, promoters_data = ensembl_promoters)
mmp13_loops <- get_gene_loops(gene = "MMP13", loops = CQTL_10kb_sig_gained_loops, 
                              loop_binsize = 30000, promoters_data = ensembl_promoters)
# DATA FOR FNF RESPONSE EGENE WTIH LOOP -------------------------------

lpar1_loop_egene <- read_csv("/proj/phanstiel_lab/Data/processed/CQTL/eqtl/FNF_PEER_k22_genoPC_cond1Mb_topSignals_rsID_signalRanges.csv") |> 
  filter(gene_symbol == "LPAR1")

lpar1_loop_egene_ld <- fread("/proj/phanstiel_lab/Data/processed/CQTL/eqtl/FNF_PEER_k22_genoPC_cond1Mb_topSignals_rsID_LD.csv",
                             data.table = FALSE) |>
  filter(gene_symbol == "LPAR1") |>
  dplyr::select(ld_variantID, R2)

lpar1_signals <- bind_rows(fread(paste0("/proj/phanstiel_lab/Data/processed/CQTL/eqtl/CTL_PEER_k20_genoPC_allSignals_nom1Mb_MAFs_",
                                        lpar1_loop_egene[["gene_chr"]], ".csv"),
                                  data.table = FALSE) |>
                              filter(gene_id == lpar1_loop_egene[["gene_id"]]) |>
                              mutate(Condition = "PBS"),
                            fread(paste0("/proj/phanstiel_lab/Data/processed/CQTL/eqtl/FNF_PEER_k22_genoPC_allSignals_nom1Mb_MAFs_",
                                         lpar1_loop_egene[["gene_chr"]], ".csv"),
                                  data.table = FALSE) |>
                              filter(gene_id == lpar1_loop_egene[["gene_id"]]) |>
                              mutate(Condition = "FN-f")) |>
  mutate(rsID = ifelse(variantID == lpar1_loop_egene[["variantID"]], lpar1_loop_egene[["rsID"]], NA)) |>
  dplyr::rename(chrom = variant_chr,
                pos = variant_start,
                p = nom_pval,
                snp = rsID) |>
  left_join(lpar1_loop_egene_ld, by = join_by(variantID == ld_variantID)) |>
  mutate(LDgrp = cut(R2, c(0, 0.2, 0.4, 0.6, 0.8, 1)))
lpar1_signals$LDgrp <- addNA(lpar1_signals$LDgrp)

# Get loop

# Ensembl TxDb
ensembl_txdb <- makeTxDbFromEnsembl(organism = "Homo sapiens", release = 111)
seqlevelsStyle(ensembl_txdb) <- "UCSC"

load("/proj/phanstiel_lab/Data/processed/CQTL/hic/loops/CQTL_10kb_static_loops.rda")
ensembl_promoter_genes <- promoters(ensembl_txdb, columns = c("TXID", "TXNAME", "GENEID"))
# convert GENEID CharacterList just to character
ensembl_promoter_genes$GENEID <- as.character(ensembl_promoter_genes$GENEID)

lpar1_loop_10kb <- get_signal_eGene_loops(eSNP_eGene = lpar1_loop_egene, 
                                    promoters_mapped_genes = ensembl_promoter_genes, 
                                    loops = CQTL_10kb_static_loops,
                                    loop_binsize = 30000)
lpar1_loop_30kb <- snapToBins(lpar1_loop_10kb, binSize = 30000)
# plotgardener layout -----------------------------------------------------

pdf("plots/distal_Fig4/Fig4_v2.pdf", width = 11, height = 6.5)
pageCreate(width = 11, height = 6.5, showGuides = FALSE)

### A - Jun
plotText("A", x = 0.05, y = 0.1, just = c("left", "top"), fontfamily = "Helvetica",
         fontsize = 14, fontface = "bold")
jun_loop_region <- pgParams(assembly = "hg38",
                            chrom = seqnames1(jun_loops),
                            chromstart = start1(jun_loops) - 100000,
                            chromend = end2(jun_loops) + 100000, 
                            x = 0.25, width = 2, norm = "SCALE")

# PBS Hi-C
jun_pbs_hic <- plotHicRectangle(data = "/proj/phanstiel_lab/Data/processed/CQTL/hic/mergedAll/PBS_inter_30.hic", 
                                resolution = 10000, params = jun_loop_region,
                                y = 0.15, height = 1, zrange = c(0, 100))
annoPixels(plot = jun_pbs_hic, data = jun_loops, type = "arrow", shift = 6)
annoHeatmapLegend(plot = jun_pbs_hic, width = 0.065, height = 0.6, fontsize = 5,
                  x = 2.3, y = 0.15)
plotText(label = "PBS", fontsize = 7, fontfamily = "Helvetica",
         x = 0.375, y = 0.23)
# FN-f Hi-C
jun_fnf_hic <- plotHicRectangle(data = "/proj/phanstiel_lab/Data/processed/CQTL/hic/mergedAll/FNF_inter_30.hic",
                                resolution = 10000, params = jun_loop_region,
                                y = 1.2, height = 1, zrange = c(0, 100))
annoPixels(plot = jun_fnf_hic, data = jun_loops, type = "arrow", shift = 6)
plotText(label = "FN-f", fontsize = 7, fontfamily = "Helvetica",
         x = 0.375, y = 1.28)

# ATAC
plotText(label = "ATAC", rot = 90, fontcolor = yl_gn_bu[4],
         fontsize = 7, x = 0.15, y = 2.5, fontfamily = "Helvetica")
plotMultiSignal(data = list("/proj/phanstiel_lab/Data/processed/CQTL/atac/atac_combined_n3/mergedSignal/CTL_merged.bw",
                            "/proj/phanstiel_lab/Data/processed/CQTL/atac/atac_combined_n3/mergedSignal/FNF_merged.bw"),
                params = jun_loop_region, 
                y = 2.25, height = 0.5, linecolor = yl_gn_bu[4], 
                fill = yl_gn_bu[4], gapdistance = 0.1)
plotText(label = "PBS", fontcolor = yl_gn_bu[4],
         fontsize = 6, x = 0.25, y = 2.275, just = "left", fontfamily = "Helvetica")
plotText(label = "FN-f", fontcolor = yl_gn_bu[4],
         fontsize = 6, x = 0.25, y = 2.5, just = "left", fontfamily = "Helvetica")

# RNA
plotText(label = "RNA", rot = 90, fontcolor = yl_gn_bu[8],
         fontsize = 7, x = 0.15, y = 3.05)
plotMultiSignal(data = list("/proj/phanstiel_lab/Data/processed/CQTL/rna/CQTL_AM7180_AM7634/mergeSignal/CQTL_CTL.bw",
                            "/proj/phanstiel_lab/Data/processed/CQTL/rna/CQTL_AM7180_AM7634/mergeSignal/CQTL_FNF.bw"),
                params = jun_loop_region, 
                y = 2.775, height = 0.5, linecolor = yl_gn_bu[8], 
                fill = yl_gn_bu[8], gapdistance = 0.05)
plotText(label = "PBS", fontcolor = yl_gn_bu[8],
         fontsize = 6, x = 0.25, y = 2.825, just = "left")
plotText(label = "FN-f", fontcolor = yl_gn_bu[8],
         fontsize = 6, x = 0.25, y = 3.075, just = "left")

# Genes
jun_gene_plot <- plotGenes(params = jun_loop_region, 
                             geneHighlights = data.frame("gene" = "JUN",
                                                         "color" = yl_gn_bu[6]),
                             y = 3.3, height = 0.4, fontsize = 8,
                           stroke = 0.05)
annoGenomeLabel(plot = jun_gene_plot, y = 3.725, params = jun_loop_region,
                fontsize = 6, fontcolor = "grey20", linecolor = "grey20")

### B - IL6
plotText("B", x = 2.475, y = 0.1, just = c("left", "top"), fontfamily = "Helvetica",
         fontsize = 14, fontface = "bold")

il6_loop_region <- pgParams(assembly = "hg38",
                            chrom = unique(seqnames1(il6_loops)),
                            chromstart = min(start1(il6_loops)) - 100000,
                            chromend = max(end2(il6_loops)) + 100000, 
                            x = 2.65, width = 2, norm = "SCALE")

# PBS Hi-C
il6_pbs_hic <- plotHicRectangle(data = "/proj/phanstiel_lab/Data/processed/CQTL/hic/mergedAll/PBS_inter_30.hic", 
                                resolution = 10000, params = il6_loop_region,
                                y = 0.15, height = 1, zrange = c(0, 140))
annoPixels(plot = il6_pbs_hic, data = il6_loops[1], type = "arrow", shift = 4)
annoHeatmapLegend(plot = il6_pbs_hic, width = 0.065, height = 0.6, fontsize = 5,
                  x = 4.7, y = 0.15)
plotText(label = "PBS", fontsize = 7, fontfamily = "Helvetica",
         x = 2.775, y = 0.23)

# FN-f Hi-C
il6_fnf_hic <- plotHicRectangle(data = "/proj/phanstiel_lab/Data/processed/CQTL/hic/mergedAll/FNF_inter_30.hic",
                                resolution = 10000, params = il6_loop_region,
                                y = 1.2, height = 1, zrange = c(0, 140))
annoPixels(plot = il6_fnf_hic, data =  il6_loops[1], type = "arrow", shift = 4)
plotText(label = "FN-f", fontsize = 7, fontfamily = "Helvetica",
         x = 2.775, y = 1.28)

# ATAC
plotMultiSignal(data = list("/proj/phanstiel_lab/Data/processed/CQTL/atac/atac_combined_n3/mergedSignal/CTL_merged.bw",
                            "/proj/phanstiel_lab/Data/processed/CQTL/atac/atac_combined_n3/mergedSignal/FNF_merged.bw"),
                params = il6_loop_region, 
                y = 2.25, height = 0.5, linecolor = yl_gn_bu[4], 
                fill = yl_gn_bu[4], gapdistance = 0.1)
plotText(label = "PBS", fontcolor = yl_gn_bu[4],
         fontsize = 6, x = 2.65, y = 2.275, just = "left", fontfamily = "Helvetica")
plotText(label = "FN-f", fontcolor = yl_gn_bu[4],
         fontsize = 6, x = 2.65, y = 2.5, just = "left", fontfamily = "Helvetica")

# RNA
plotMultiSignal(data = list("/proj/phanstiel_lab/Data/processed/CQTL/rna/CQTL_AM7180_AM7634/mergeSignal/CQTL_CTL.bw",
                            "/proj/phanstiel_lab/Data/processed/CQTL/rna/CQTL_AM7180_AM7634/mergeSignal/CQTL_FNF.bw"),
                params = il6_loop_region, 
                y = 2.775, height = 0.5, linecolor = yl_gn_bu[8], 
                fill = yl_gn_bu[8], gapdistance = 0.05)
plotText(label = "PBS", fontcolor = yl_gn_bu[8],
         fontsize = 6, x = 2.65, y = 2.825, just = "left")
plotText(label = "FN-f", fontcolor = yl_gn_bu[8],
         fontsize = 6, x = 2.65, y = 3.075, just = "left")

# Genes
il6_gene_plot <- plotGenes(params = il6_loop_region, 
                           geneHighlights = data.frame("gene" = "IL6",
                                                       "color" = yl_gn_bu[6]),
                           y = 3.3, height = 0.4, fontsize = 8,
                           stroke = 0.05, strandLabels = FALSE)
annoGenomeLabel(plot = il6_gene_plot, y = 3.725, params = il6_loop_region,
                fontsize = 6, fontcolor = "grey20", linecolor = "grey20")

## C - MMP13
plotText("C", x = 4.9, y = 0.1, just = c("left", "top"), fontfamily = "Helvetica",
         fontsize = 14, fontface = "bold")

mmp13_loop_region <- pgParams(assembly = "hg38",
                            chrom = seqnames1(mmp13_loops),
                            chromstart = start1(mmp13_loops) - 90000,
                            chromend = end2(mmp13_loops) + 110000, 
                            x = 5.1, width = 2, norm = "SCALE")
# PBS Hi-C
mmp13_pbs_hic <- plotHicRectangle(data = "/proj/phanstiel_lab/Data/processed/CQTL/hic/mergedAll/PBS_inter_30.hic", 
                                resolution = 10000, params = mmp13_loop_region,
                                y = 0.15, height = 1, zrange = c(0, 150))
annoPixels(plot = mmp13_pbs_hic, data = mmp13_loops, type = "arrow", shift = 4)
annoHeatmapLegend(plot = mmp13_pbs_hic, width = 0.065, height = 0.6, fontsize = 5,
                  x = 7.15, y = 0.15)
plotText(label = "PBS", fontsize = 7, fontfamily = "Helvetica",
         x = 5.225, y = 0.23)
# FN-f Hi-C
mmp13_fnf_hic <- plotHicRectangle(data = "/proj/phanstiel_lab/Data/processed/CQTL/hic/mergedAll/FNF_inter_30.hic",
                                resolution = 10000, params = mmp13_loop_region,
                                y = 1.2, height = 1, zrange = c(0, 150))
annoPixels(plot = mmp13_fnf_hic, data = mmp13_loops, type = "arrow", shift = 4)
plotText(label = "FN-f", fontsize = 7, fontfamily = "Helvetica",
         x = 5.225, y = 1.28)

# ATAC
plotMultiSignal(data = list("/proj/phanstiel_lab/Data/processed/CQTL/atac/atac_combined_n3/mergedSignal/CTL_merged.bw",
                            "/proj/phanstiel_lab/Data/processed/CQTL/atac/atac_combined_n3/mergedSignal/FNF_merged.bw"),
                params = mmp13_loop_region, 
                y = 2.25, height = 0.5, linecolor = yl_gn_bu[4], 
                fill = yl_gn_bu[4], gapdistance = 0.1)
plotText(label = "PBS", fontcolor = yl_gn_bu[4],
         fontsize = 6, x = 5.1, y = 2.275, just = "left", fontfamily = "Helvetica")
plotText(label = "FN-f", fontcolor = yl_gn_bu[4],
         fontsize = 6, x = 5.1, y = 2.5, just = "left", fontfamily = "Helvetica")

# RNA
plotMultiSignal(data = list("/proj/phanstiel_lab/Data/processed/CQTL/rna/CQTL_AM7180_AM7634/mergeSignal/CQTL_CTL.bw",
                            "/proj/phanstiel_lab/Data/processed/CQTL/rna/CQTL_AM7180_AM7634/mergeSignal/CQTL_FNF.bw"),
                params = mmp13_loop_region, 
                y = 2.775, height = 0.5, linecolor = yl_gn_bu[8], 
                fill = yl_gn_bu[8], gapdistance = 0.05)
plotText(label = "PBS", fontcolor = yl_gn_bu[8],
         fontsize = 6, x = 5.1, y = 2.825, just = "left")
plotText(label = "FN-f", fontcolor = yl_gn_bu[8],
         fontsize = 6, x = 5.1, y = 3.075, just = "left")

# Genes
mmp13_gene_plot <- plotGenes(params = mmp13_loop_region, 
                           geneHighlights = data.frame("gene" = "MMP13",
                                                       "color" = yl_gn_bu[6]),
                           y = 3.3, height = 0.4, fontsize = 8,
                           stroke = 0.05, strandLabels = FALSE)
annoGenomeLabel(plot = mmp13_gene_plot, y = 3.725, params = mmp13_loop_region,
                fontsize = 6, fontcolor = "grey20", linecolor = "grey20")


## D - density of genes at static and gained loop anchors
plotText("D", x = 7.35, y = 0.1, just = c("left", "top"), fontfamily = "Helvetica",
         fontsize = 14, fontface = "bold")
plotGG(plot = anchor_expression_density, 
       x = 7.5, y = 0, width = 3.5, height = 2)

## E - barplot of GO enrichment terms for gained anchor genes
plotText("E", x = 0.05, y = 4, just = c("left", "top"), fontfamily = "Helvetica",
         fontsize = 14, fontface = "bold")
plotGG(plot = loop_go_barplot, x = 0.2, y = 3.95, width = 2.85, height = 2.55)

## F - distal/loop barplot 
plotText("F", x = 3.05, y = 4, just = c("left", "top"), fontfamily = "Helvetica",
         fontsize = 14, fontface = "bold")
plotGG(plot = distal_ld_pixel_contactfreq_boxplot, 
       x = 3.1, y = 4, 
       width = 2.1, height = 2.5)

## G - distal lead SNP contact freq boxplot
plotText("G", x = 5, y = 4, just = c("left", "top"), fontfamily = "Helvetica",
         fontsize = 14, fontface = "bold")
plotGG(eqtl_cat_distal_contactfreq_boxplot,
       x = 5.1, y = 4, width = 2.15, height = 2.45)

## H - example region
plotText("H", x = 7.35, y = 1.95, just = c("left", "bottom"), 
         fontfamily = "Helvetica",
         fontsize = 14, fontface = "bold")

lpar1_region <- pgParams(assembly = "hg38",
                        chrom = lpar1_loop_egene[["gene_chr"]],
                           chromstart = as.numeric(lpar1_loop_egene[["variant_start"]]) - 400000,
                           chromend = as.numeric(lpar1_loop_egene[["variant_start"]]) + 500000,
                          x = 7.7, width = 3.1)

lpar1_ylim <- ceiling(max(log10(lpar1_signals$p)*-1)) + 2

pbs_lpar1_threshold <- lpar1_signals |> 
  filter(Condition == "PBS") |> 
  pull(pval_nominal_threshold) |> 
  unique()

fnf_lpar1_threshold <- lpar1_signals |> 
  filter(Condition == "FN-f") |> 
  pull(pval_nominal_threshold) |> 
  unique()

# PBS eQTL
pbs_lpar1 <- plotManhattan(data = lpar1_signals |> filter(Condition == "PBS"),
                            params = lpar1_region, 
                            range = c(0, lpar1_ylim),
                            snpHighlights = data.frame(snp = lpar1_loop_egene[["rsID"]],
                                                       pch = 23,
                                                       cex = 0.5),
                            fill = colorby("LDgrp",
                                           palette = colorRampPalette(c("#262C74",
                                                                        "#98CDED",
                                                                        "#499A53",
                                                                        "#EEA741",
                                                                        "#DD3931",
                                                                        "grey"))),
                            y = 2, height = 1)
annoSegments(plot = pbs_lpar1, x0 = unit(0, "npc"), x1 = unit(1, "npc"),
             y0 = -1*log10(pbs_lpar1_threshold), y1 = -1*log10(pbs_lpar1_threshold),
             default.units = "native", lty = 2)

annoYaxis(plot = pbs_lpar1, at = seq(0, lpar1_ylim, 2), 
          axisLine = TRUE, fontsize = 6)
plotText(
  label = "-log10(p-value)", x = 7.425, y = 2.5, rot = 90,
  fontsize = 6, just = "center",
  default.units = "inches", fontfamily = "Helvetica"
)
plotText(label = "PBS", x = 7.725, y = 2, just = c("left", "top"),
         fontfamily = "Helvetica", fontsize = 10)

plotLegend(legend = c("0.8 - 1.0",
                      "0.6 - 0.8",
                      "0.4 - 0.6",
                      "0.2 - 0.4",
                      "0.0 - 0.2"),
           fill = c("#DD3931", "#EEA741", "#499A53","#98CDED","#262C74"),
           x = 10, y = 2.05, width = 0.1, height = 0.35, border = FALSE, 
           fontsize = 6)
plotText("r2", x = 10.55, y = 1.975, fontsize = 6, fontfamily = "Helvetica")

# FN-f eQTL
fnf_lpar1 <- plotManhattan(data = lpar1_signals |> filter(Condition == "FN-f"),
                              params = lpar1_region, 
                              range = c(0, lpar1_ylim),
                              snpHighlights = data.frame(snp = lpar1_loop_egene[["rsID"]],
                                                         pch = 23,
                                                         cex = 0.5),
                              fill = colorby("LDgrp",
                                             palette = colorRampPalette(c("#262C74",
                                                                          "#98CDED",
                                                                          "#499A53",
                                                                          "#EEA741",
                                                                          "#DD3931",
                                                                          "grey"))),
                              y = 3.1, height = 1)
annoSegments(plot = fnf_lpar1, x0 = unit(0, "npc"), x1 = unit(1, "npc"),
             y0 = -1*log10(fnf_lpar1_threshold), y1 = -1*log10(fnf_lpar1_threshold),
             default.units = "native", lty = 2)
annoYaxis(plot = fnf_lpar1, at = seq(0, lpar1_ylim, 2), 
          axisLine = TRUE, fontsize = 6)
plotText(
  label = "-log10(p-value)", x = 7.425, y = 3.6, rot = 90,
  fontsize = 6, just = "center",
  default.units = "inches", fontfamily = "Helvetica"
)

plotText(label = "FN-f", x = 7.725, y = 3.1, just = c("left", "top"),
         fontfamily = "Helvetica", fontsize = 10)

plotText(label = lpar1_loop_egene[["rsID"]], x = 9.3, y = 3.175,
         fontfamily = "Helvetica", fontsize = 7)


# Hi-C with megamap for static loop

lpar1_hic <- plotHicRectangle(data = "/proj/phanstiel_lab/Data/processed/CQTL/hic/mergedAll/CQTL_inter_30.hic",
                              params = lpar1_region,
                              norm = "SCALE", resolution = 10000,
                              y = 4.15, height = 1,
                              zrange = c(0, 175))

annoHeatmapLegend(
  plot = lpar1_hic,
  x = 10.825,
  y = 4.15,
  width = 0.065, height = 0.6,
  fontsize = 5
)

annoPixels(plot = lpar1_hic, data = lpar1_loop_10kb, type = "arrow",
           shift = 6)

# ATAC
plotText(label = "ATAC", rot = 90, fontcolor = yl_gn_bu[4],
         fontsize = 7, x = 7.6, y = 5.3, fontfamily = "Helvetica")
plotMultiSignal(data = list("/proj/phanstiel_lab/Data/processed/CQTL/atac/atac_combined_n3/mergedSignal/CTL_merged.bw",
                            "/proj/phanstiel_lab/Data/processed/CQTL/atac/atac_combined_n3/mergedSignal/FNF_merged.bw"),
                params = lpar1_region, 
                y = 5.2, height = 0.3, linecolor = yl_gn_bu[4], 
                fill = yl_gn_bu[4], gapdistance = 0.05)
plotText(label = "PBS", fontcolor = yl_gn_bu[4],
         fontsize = 6, x = 7.7, y = 5.2, just = "left", fontfamily = "Helvetica")
plotText(label = "FN-f", fontcolor = yl_gn_bu[4],
         fontsize = 6, x = 7.7, y = 5.375, just = "left", fontfamily = "Helvetica")

# RNA
plotText(label = "RNA", rot = 90, fontcolor = yl_gn_bu[8],
         fontsize = 7, x = 7.6, y = 5.7)
plotMultiSignal(data = list("/proj/phanstiel_lab/Data/processed/CQTL/rna/CQTL_AM7180_AM7634/mergeSignal/CQTL_CTL.bw",
                            "/proj/phanstiel_lab/Data/processed/CQTL/rna/CQTL_AM7180_AM7634/mergeSignal/CQTL_FNF.bw"),
                params = lpar1_region, 
                y = 5.55, height = 0.3, linecolor = yl_gn_bu[8], 
                fill = yl_gn_bu[8], gapdistance = 0.05)
plotText(label = "PBS", fontcolor = yl_gn_bu[8],
         fontsize = 6, x = 7.7, y = 5.55, just = "left")
plotText(label = "FN-f", fontcolor = yl_gn_bu[8],
         fontsize = 6, x = 7.7, y = 5.725, just = "left")

# Genes
lpar1_genes <- plotGenes(params = lpar1_region, y = 5.85,
                         height = 0.4, geneHighlights = data.frame("gene" = "LPAR1",
                                                                   "color" = "#37a7db"),
                        fontsize = 6)
annoGenomeLabel(plot = lpar1_genes, y = 6.25, params = lpar1_region,
                fontsize = 6, fontcolor = "grey20", linecolor = "grey20")

dev.off()
