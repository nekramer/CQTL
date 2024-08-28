library(tidyverse)
library(plotgardener)
source("../plotting_utils.R")

# Functions ---------------------------------------------------------------

get_skipped_transcripts_locus <- function(lead_eSNP_eGene, ld_data, gene_promoters){
  
  eGene_id <- lead_eSNP_eGene[["gene_id"]]
  print(eGene_id)
  
  # Grab LD buddies from ld_data
  ld_data_locus <- ld_data |> 
    filter(rsID == lead_eSNP_eGene[["rsID"]] &
             gene_id == lead_eSNP_eGene[["gene_id"]] &
             signal == lead_eSNP_eGene[["signal"]]) |> 
    separate_wider_delim(cols = "ld_variantID",
                         delim = ":",
                         names = c(NA, "ld_pos", NA, NA), cols_remove = FALSE)
  
  # Define locus range from LD buddies
  locus_range <- GRanges(seqnames = unique(ld_data_locus$gene_chr),
                         ranges = IRanges(start = min(as.numeric(ld_data_locus$ld_pos)),
                                          end = max(as.numeric(ld_data_locus$ld_pos))))
  
  # Get distance between locus and eGene
  egene_distance <- distance(locus_range, GRanges(seqnames = unique(ld_data_locus$gene_chr),
                                                  ranges = IRanges(start = unique(ld_data_locus$gene_start),
                                                                   end = unique(ld_data_locus$gene_end)),
                                                  strand = unique(ld_data_locus$gene_strand)))
  
  
  # Grab genes with a promoter within 1 Mb of the locus
  genes_1mb_hits <- findOverlaps(locus_range, gene_promoters, maxgap = 1e6)
  # Get genes with distance within 1 Mb
  genes_1mb <- gene_promoters[unique(subjectHits(genes_1mb_hits))]
  # Calculate distances 
  genes_1mb$locus_dist <- distance(locus_range, genes_1mb) 
  
  
  # Order by distance
  gene_1mb_ordered <- genes_1mb |> 
    as_tibble() |> 
    dplyr::select(GENEID, locus_dist) |> 
    bind_rows(tibble(GENEID = eGene_id, locus_dist = egene_distance)) |> 
    distinct() |> 
    arrange(locus_dist)
  
  # Calculate number of genes skipped by finding eGene id in ordered list, grab first one
  egene_index <- which(gene_1mb_ordered$GENEID == eGene_id)[1]
  
  # Grab all the genes in the list before that
  locus_skipped_genes <- gene_1mb_ordered[1:egene_index-1,] |> 
      pull(GENEID) |> 
      unique() |> 
      length()
  
  return(locus_skipped_genes)
  
}


get_skipped_transcripts_lead <- function(lead_eSNP_eGene, ld_data, gene_promoters){
  
  eGene_id <- lead_eSNP_eGene[["gene_id"]]
  print(eGene_id)
  
  
  # Define locus range from lead
  locus_range <- GRanges(seqnames = lead_eSNP_eGene[["gene_chr"]],
                         ranges = IRanges(start = as.numeric(lead_eSNP_eGene[["variant_start"]]),
                                          end = as.numeric(lead_eSNP_eGene[["variant_start"]])))
  
  # Get distance between locus and eGene
  egene_distance <- distance(locus_range, GRanges(seqnames = lead_eSNP_eGene[["gene_chr"]],
                                                  ranges = IRanges(start = as.numeric(lead_eSNP_eGene[["gene_start"]]),
                                                                   end = as.numeric(lead_eSNP_eGene[["gene_end"]])),
                                                  strand = lead_eSNP_eGene[["gene_strand"]]))
  
  
  # Grab genes with a promoter within 1 Mb of the locus
  genes_1mb_hits <- findOverlaps(locus_range, gene_promoters, maxgap = 1e6)
  # Get genes with distance within 1 Mb
  genes_1mb <- gene_promoters[unique(subjectHits(genes_1mb_hits))]
  # Calculate distances 
  genes_1mb$locus_dist <- distance(locus_range, genes_1mb) 
  
  
  # Order by distance
  gene_1mb_ordered <- genes_1mb |> 
    as_tibble() |> 
    dplyr::select(GENEID, locus_dist) |> 
    bind_rows(tibble(GENEID = eGene_id, locus_dist = egene_distance)) |> 
    distinct() |> 
    arrange(locus_dist)
  
  # Calculate number of genes skipped by finding eGene id in ordered list, grab first one
  egene_index <- which(gene_1mb_ordered$GENEID == eGene_id)[1]
  
  # Grab all the genes in the list before that
  locus_skipped_genes <- gene_1mb_ordered[1:egene_index-1,] |> 
    pull(GENEID) |> 
    unique() |> 
    length()
  
  return(locus_skipped_genes)
  
}
# Function to get the transcripts with promoters that overlap an eSNP_eGene signal range
# @returns transcript IDs
# get_skipped_transcripts <- function(eSNP_eGene, gene_promoters){
#   snp_pos <- as.numeric(eSNP_eGene[["variant_start"]])
#   gene_tss <- as.numeric(eSNP_eGene[["gene_start"]])
#   strand <- eSNP_eGene[["gene_strand"]]
#   chrom <- eSNP_eGene[["gene_chr"]]
#   
#   snp_gene_range <- GRanges(seqnames = chrom, 
#                             ranges = IRanges(start = min(c(snp_pos, gene_tss)), 
#                                              end = max(c(snp_pos, gene_tss))),
#                             strand = strand)
#   overlappingPromoters <- subsetByOverlaps(gene_promoters, snp_gene_range)
#   
#   return(overlappingPromoters$tx_id)
# }


# Function to determine which eSNP_eGene signals are supported by an ABC connection
# @returns a character, either "abc" or "no_abc"
determine_abc_support <- function(eSNP_eGene, abc){
  
  # Get gene_symbol of eGene to look up in ABC data
  eGene_symbol <- eSNP_eGene[["gene_symbol"]]
  
  if (!is.na(eGene_symbol)){
    # Extract signal chrom, start, end, and strand
    signal_start <- as.numeric(eSNP_eGene[["min_signal_pos"]])
    signal_end <- as.numeric(eSNP_eGene[["max_signal_pos"]])
    strand <- eSNP_eGene[["gene_strand"]]
    chrom <- eSNP_eGene[["gene_chr"]]
    
    # Define signal range in a GRanges object
    signal_range <- GRanges(seqnames = chrom,
                            ranges = IRanges(start = signal_start,
                                             end = signal_end),
                            strand = strand)
    
    # Filter for eGene_symbol in ABC data (only has gene symbols)
    eGene_abc <- abc |> 
      filter(TargetGene == eGene_symbol)
    
    # Convert gene ABC data into GRanges
    eGene_abc_gr <- makeGRangesFromDataFrame(eGene_abc, keep.extra.columns = TRUE)
    
    # Find overlaps between signal and ABC 
    abc_links <- countOverlaps(signal_range, eGene_abc_gr)
    
    if (abc_links > 0){
      return("abc")
    } else {
      return("no_abc")
    }
    
  } else {
    return("no_abc")
  }
}

# Distribution of distance between lead SNP and eGene for each independent signal --------------------

pbs_egenes_signals <- read_csv("/proj/phanstiel_lab/Data/processed/CQTL/eqtl/CTL_PEER_k20_genoPC_cond1Mb_topSignals_rsID.csv") |> 
  mutate(tss_dist = abs(gene_start - variant_start))

fnf_egenes_signals <- read_csv("/proj/phanstiel_lab/Data/processed/CQTL/eqtl/FNF_PEER_k22_genoPC_cond1Mb_topSignals_rsID.csv") |> 
  mutate(tss_dist = abs(gene_start - variant_start))

all_egenes_signals <- bind_rows(pbs_egenes_signals, fnf_egenes_signals) |> 
  distinct(rsID, gene_id, .keep_all = TRUE)

# Percentage more than 100 kb
perc_100kb <- (all_egenes_signals |> filter(tss_dist > 100000) |> nrow())/nrow(all_egenes_signals)*100


tss_dist_plot <- ggplot(all_egenes_signals, aes(x = tss_dist)) +
  geom_line(stat = "bin", binwidth = 50000) +
  scale_x_continuous(limits = c(0, 1e6), name = "Distance to eGene TSS (Kb)", expand = c(0,0),
                     labels = label_number(scale = .001, big.mark = "")) +
  scale_y_continuous(name = "Number of lead eSNPs", expand = c(0,0),
                     limits = c(0, 2500)) +
  theme(text = element_text(family = "Helvetica"),
        panel.background = element_rect(fill = 'transparent', color = "transparent"),
        plot.background = element_rect(fill = 'transparent', color = "transparent"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(), 
        axis.ticks.x = element_blank(),
        axis.title = element_text(size = 9),
        axis.line = element_line(color = "black", 
                                 linewidth = 0.25),
        axis.ticks.length.y = unit(-0.1, "cm"),
        axis.ticks.y = element_line(color = "black",
                                    linewidth = 0.25),
        axis.text = element_text(color = "black", size = 8),
        plot.margin = margin(t = 5, r = 10, b = 5, l = 5))

y_100kb <- ggplot_build(tss_dist_plot)$data[[1]] |> 
  filter(x == 100000) |> 
  pull(count)


tss_dist_plot <- tss_dist_plot +
  annotate(geom = "point", x = 100000, y = y_100kb, size = 2.5) +
  annotate(geom = "text", label = paste0(signif(perc_100kb, digits = 3), 
                                         "% of lead eQTLs are more\nthan 100 Kb from their eGene TSS"), 
           x = 125000, y = 750, hjust = 0, size = 3.5)

save(tss_dist_plot, file = "plots/distal_SupFig6/tss_dist_plot.rda")


# Do high conf condition specific leads tend to be farther from their eGene TSS? --------

pbs_highconf <- read_csv("data/reqtl/CTL_sig01_beta_donor_reQTLs_PEER_k20_genoPC.csv")
fnf_highconf <- read_csv("data/reqtl/FNF_sig01_beta_donor_reQTLs_PEER_k22_genoPC.csv")


pbs_egenes_signals <- read_csv("data/eqtl/CTL_PEER_k20_genoPC_cond1Mb_topSignals_rsID.csv") |> 
  mutate(tss_dist = abs(gene_start - variant_start),
         response = ifelse(gene_id %in% pbs_highconf$gene_id, "condition-specific", "not condition-specific"))

fnf_egenes_signals <- read_csv("data/eqtl/FNF_PEER_k22_genoPC_cond1Mb_topSignals_rsID.csv") |> 
  mutate(tss_dist = abs(gene_start - variant_start),
         response = ifelse(gene_id %in% fnf_highconf$gene_id,  "condition-specific", "not condition-specific"))

all_egenes_signals <- bind_rows(pbs_egenes_signals, fnf_egenes_signals) |> 
  distinct(rsID, gene_id, .keep_all = TRUE)


wilcox.test(x = all_egenes_signals |> filter(response == "condition-specific") |> pull(tss_dist),
            y = all_egenes_signals |> filter(response == "not condition-specific") |> pull(tss_dist),
            alternative = "greater")


ggplot(all_egenes_signals, mapping = aes(x = response, y = tss_dist)) +
  geom_boxplot(aes(color = response, fill = response), outlier.shape = NA) +
  scale_y_continuous(name = "Distance to eGene TSS (Kb)", 
                     labels = label_number(scale = .001, big.mark = ""),
                     limits = c(0, 500000)) +
  scale_fill_manual(values = c("#FBBE67", "grey80")) + 
  scale_color_manual(values = c("#97723e", "grey20")) + 
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.y = element_line(color = "black", linewidth = 0.2),
        axis.ticks.length.y = unit(-0.1, "cm"),
        axis.ticks.y = element_line(color = "black", linewidth = 0.2),
        panel.background = element_rect(fill = 'transparent', color = "transparent"),
        plot.background = element_rect(fill = 'transparent', color = "transparent"),
        axis.text = element_text(color = "black"),
        axis.text.x = element_markdown(margin = margin(t = 3), family = "Helvetica", size = 6),
        text = element_text(family = "Helvetica", size = 8),
        axis.title.y = element_markdown(size = 7, margin = margin(r = -2))) +
  coord_cartesian(clip = "off")

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
                                                       0.2, 0.4, 0.6, 0.8))
  
# Number of genes skipped -------------------------------------------------

pbs_leads <- read_csv("/proj/phanstiel_lab/Data/processed/CQTL/eqtl/CTL_PEER_k20_genoPC_cond1Mb_topSignals_rsID.csv")
pbs_ld08 <- fread("/proj/phanstiel_lab/Data/processed/CQTL/eqtl/CTL_PEER_k20_genoPC_cond1Mb_topSignals_rsID_LD.csv",
                     data.table = FALSE) |> 
  filter(R2 >= 0.8)

  
fnf_leads <- read_csv("/proj/phanstiel_lab/Data/processed/CQTL/eqtl/FNF_PEER_k22_genoPC_cond1Mb_topSignals_rsID.csv")
fnf_ld08 <- fread("/proj/phanstiel_lab/Data/processed/CQTL/eqtl/FNF_PEER_k22_genoPC_cond1Mb_topSignals_rsID_LD.csv",
                  data.table = FALSE) |> 
  filter(R2 >= 0.8)


# Get ENSEMBL Txdb
ensembl_txdb <- makeTxDbFromEnsembl(organism = "Homo sapiens", release = 111)
# Update seqlevels to UCSC style for compatibility with eqtl data
seqlevelsStyle(ensembl_txdb) <- "UCSC"

ensembl_promoters <- promoters(ensembl_txdb)

# Map transcripts IDs to gene IDs with TxDb
tx_genes <- AnnotationDbi::select(ensembl_txdb, 
                                  keys = as.character(ensembl_promoters$tx_id),
                                  keytype = "TXID",
                                  columns = "GENEID")

ensembl_promoters_genes <- ensembl_promoters |> 
  as_tibble() |> 
  left_join(tx_genes, by = join_by(tx_id == TXID)) |> 
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)
  
# Get number of skipped genes per signal locus with LD > 0.8

## PBS

pbs_skipped_locus <- apply(pbs_leads, 1, get_skipped_transcripts_locus, ld_data = pbs_ld08,
      gene_promoters = ensembl_promoters_genes)

## FN-f

fnf_skipped_locus <- apply(fnf_leads, 1, get_skipped_transcripts_locus, ld_data = fnf_ld08,
                           gene_promoters = ensembl_promoters_genes)

# Combine for plotting

all_skipped_locus <- bind_rows(tibble(Condition = "PBS",
                                      rsID = pbs_leads$rsID,
                                      gene_id = pbs_leads$gene_id,
                                      gene_symbol = pbs_leads$gene_symbol,
                                      signal = pbs_leads$signal,
                                      skipped_genes = pbs_skipped_locus),
                               tibble(Condition = "FN-f",
                                      rsID = fnf_leads$rsID,
                                      gene_id = fnf_leads$gene_id,
                                      gene_symbol = fnf_leads$gene_symbol,
                                      signal = fnf_leads$signal,
                                      skipped_genes = fnf_skipped_locus)) |> 
  mutate(skipped_gene_group = ifelse(skipped_genes >= 10, "10+", as.character(skipped_genes))) |> 
  mutate(skipped_gene_group = factor(skipped_gene_group,
                                     levels = c(as.character(seq(0,9)), "10+")))


percent_1plus_skipped_genes <- (all_skipped_locus |> filter(skipped_genes >= 1) |> nrow())/nrow(all_skipped_locus) * 100 


skipped_locus_barplot <- ggplot(all_skipped_locus, aes(x = skipped_gene_group)) +
  geom_bar(stat = "count") +
  annotate(geom = "segment", x = "1", xend = "10+", y = 1250, yend = 1250,
           linewidth = 0.4) + 
  annotate(geom = "segment", x = "1", xend = "1", y = 1250, yend = 1150,
           linewidth = 0.4) +
  annotate(geom = "segment", x = "10+", xend = "10+", y = 1250, yend = 1150,
           linewidth = 0.4) +
  annotate(geom = "text", x = 6.5, y = 1650, label = paste0(signif(percent_1plus_skipped_genes, digits = 3), 
                                                            "% of eQTL loci\nskip at least one gene"),
           size = 3.5) +
  scale_x_discrete(name = "Number of genes skipped") +
  scale_y_continuous(name = "Number of eQTL loci", expand = c(0, 0)) +
  theme(text = element_text(family = "Helvetica"),
        panel.background = element_rect(fill = 'transparent', color = 'transparent'),
        plot.background = element_rect(fill = 'transparent', color = 'transparent'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        axis.ticks.x = element_blank(),
        axis.title = element_text(size = 9),
        axis.line = element_line(color = "black", 
                                 linewidth = 0.25),
        axis.ticks.length.y = unit(-0.1, "cm"),
        axis.ticks.y = element_line(color = "black",
                                    linewidth = 0.25),
        axis.text = element_text(color = "black", size = 8))

save(skipped_locus_barplot, file = "plots/distal_SupFig7/skipped_locus_barplot.rda")

# # Get number of skipped transcripts per lead eSNP-eGene pair
# pbs_signals <- read_csv("/proj/phanstiel_lab/Data/processed/CQTL/eqtl/CTL_PEER_k20_genoPC_cond1Mb_topSignals_rsID.csv")
# fnf_signals <- read_csv("/proj/phanstiel_lab/Data/processed/CQTL/eqtl/FNF_PEER_k22_genoPC_cond1Mb_topSignals_rsID.csv")
# 
# all_signals <- bind_rows(pbs_signals, fnf_signals) |>
#   distinct(rsID, gene_id, .keep_all = TRUE)
# 
# skipped_transcripts <- apply(all_signals, 1, get_skipped_transcripts_lead,
#                              gene_promoters = ensembl_promoters_genes)
# 
# skipped_transcripts <- apply(all_signals, 1, get_skipped_transcripts,
#                              gene_promoters = ensembl_promoters)
# all_signals$skipped_transcripts <- skipped_transcripts
# 
# all_signals <- all_signals |>
#   # Give each skipped transcript a separate row
#   unnest(cols = "skipped_transcripts", keep_empty = TRUE) |>
#   mutate(skipped_transcripts = as.character(skipped_transcripts))
# 
# 
# # Map transcripts IDs to gene IDs with TxDb
# tx_genes <- AnnotationDbi::select(ensembl_txdb,
#                                   keys = unique(all_signals$skipped_transcripts),
#                                   keytype = "TXID",
#                                   columns = "GENEID")
# 
# egenes_signals_skipped_genes <- left_join(all_signals, tx_genes,
#                                   by = join_by("skipped_transcripts" == "TXID"),
#                                   relationship = "many-to-many") |>
#   dplyr::select(-skipped_transcripts) |>
#   dplyr::rename(skipped_geneid = GENEID) |>
#   # Grab unique skipped genes
#   distinct() |>
#   mutate(eSNP_eGene = paste0(rsID, "_", gene_id))
# 
# 
# # Pull out ones that have a gene overlapping themself and check if they have others
# self_genes <- egenes_signals_skipped_genes |>
#   filter(gene_id == skipped_geneid)
# 
# after_self_gene_removal <- egenes_signals_skipped_genes |>
#   filter(!gene_id == skipped_geneid)
# 
# self_genes_none <- egenes_signals_skipped_genes |>
#   filter(eSNP_eGene %in% self_genes$eSNP_eGene & !eSNP_eGene %in% after_self_gene_removal$eSNP_eGene)
# 
# #Pull out ones with 0 skipped genes
# no_skipped_genes <- egenes_signals_skipped_genes |>
#   anti_join(self_genes_none) |>
#   filter(is.na(skipped_geneid))
# 
# #Group and count the rest of the skipped transcripts
# number_skipped_genes <- egenes_signals_skipped_genes |>
#   # Remove any rows where the skipped transcript maps to the eGene
#   anti_join(self_genes_none) |>
#   anti_join(no_skipped_genes) |>
#   filter(gene_id != skipped_geneid) |>
#   group_by(eSNP_eGene) |>
#   summarise(skipped_genes = dplyr::n())
# 
# 
# all_number_skipped_genes <- bind_rows(self_genes_none |>
#                                         group_by(eSNP_eGene) |>
#                                         summarise(skipped_genes = 0),
#                                       no_skipped_genes |>
#                                         group_by(eSNP_eGene) |>
#                                         summarise(skipped_genes = 0),
#                                       number_skipped_genes) |>
#   mutate(skipped_gene_group = ifelse(skipped_genes >= 10, "10+", as.character(skipped_genes))) |>
#   mutate(skipped_gene_group = factor(skipped_gene_group,
#                                      levels = c(as.character(seq(0,9)), "10+")))
# 
# 
# percent_1plus_skipped_genes <- (all_number_skipped_genes |> filter(skipped_genes >= 1) |> nrow())/nrow(all_number_skipped_genes) * 100
# 
# skipped_egene_barplot <- ggplot(all_number_skipped_genes, aes(x = skipped_gene_group)) +
#   geom_bar(stat = "count") +
#   annotate(geom = "segment", x = "1", xend = "10+", y = 1250, yend = 1250,
#            linewidth = 0.4) + 
#   annotate(geom = "segment", x = "1", xend = "1", y = 1250, yend = 1100,
#            linewidth = 0.4) +
#   annotate(geom = "segment", x = "10+", xend = "10+", y = 1250, yend = 1100,
#            linewidth = 0.4) +
#   annotate(geom = "text", x = 6.5, y = 1650, label = paste0(signif(percent_1plus_skipped_genes, digits = 3), 
#                                                             "% of lead eQTLs\nskip at least one gene"),
#            size = 3.5) +
#   scale_x_discrete(name = "Number of genes skipped") +
#   scale_y_continuous(name = "Number of lead eSNPs", expand = c(0, 0)) +
#   theme(text = element_text(family = "Helvetica"),
#         panel.background = element_rect(fill = 'transparent', color = 'transparent'),
#         plot.background = element_rect(fill = 'transparent', color = 'transparent'),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.border = element_blank(), 
#         axis.ticks.x = element_blank(),
#         axis.title = element_text(size = 9),
#         axis.line = element_line(color = "black", 
#                                  linewidth = 0.25),
#         axis.ticks.length.y = unit(-0.1, "cm"),
#         axis.ticks.y = element_line(color = "black",
#                                     linewidth = 0.25),
#         axis.text = element_text(color = "black", size = 8))
# 
# save(skipped_egene_barplot, file = "plots/distal_SupFig6/skipped_egene_barplot.rda")


# Do high conf condition specific leads tend to skip more genes ? --------

pbs_highconf <- read_csv("data/reqtl/CTL_sig01_beta_donor_reQTLs_PEER_k20_genoPC.csv")
fnf_highconf <- read_csv("data/reqtl/FNF_sig01_beta_donor_reQTLs_PEER_k22_genoPC.csv")


all_number_skipped_genes_response <- all_number_skipped_genes |> 
  separate_wider_delim(cols = "eSNP_eGene",
                       delim = "_",
                       names = c("variant", "gene_id")) |> 
  mutate(response = ifelse(gene_id %in% pbs_highconf$gene_id | 
                             gene_id %in% fnf_highconf$gene_id, TRUE, FALSE))




wilcox.test(x = all_number_skipped_genes_response |> filter(response == TRUE) |> pull(skipped_genes),
            y = all_number_skipped_genes_response |> filter(response == FALSE) |> pull(skipped_genes),
            alternative = "greater")


# Classifying distal signals as being supported by ABC or not ------------------

# Read in distal signals only
distal_all_signals_signalRanges <- read_csv("tables/eGenes_signals_distal_looping.csv")

# ABC data for PBS and FNF combined
CQTL_ABC_PBS_FNF <- read_delim("data/abc/CQTL_ABCpredictions_PBS_FNF.txt")

# Add column with abc support info
distal_all_signals_signalRanges$abc <- apply(distal_all_signals_signalRanges,
                                            1,
                                            determine_abc_support,
                                            CQTL_ABC_PBS_FNF)
# Write to file
write_csv(distal_all_signals_signalRanges,
          file = "tables/eGenes_signals_distal_looping_abc.csv")

prop_abc <- (distal_all_signals_signalRanges |> 
  filter(abc == "abc") |>
  nrow())/nrow(distal_all_signals_signalRanges)


# Assemble supp fig with plotgardener ----------------------------------------------

pdf(file = "plots/distal_SupFig7/SupFig7.pdf", width = 9.25, height = 2.5)
pageCreate(width = 9.25, height = 2.5, showGuides = FALSE)

### A
plotText("A", x = 0.2, y = 0.2, just = c("left", "top"), fontfamily = "Helvetica",
         fontsize = 14, fontface = "bold")
plotGG(plot = tss_dist_plot, x = 0.25, y = 0.25, width = 4.25, height = 2.25)


## B
plotText("B", x = 4.6, y = 0.2, just = c("left", "top"), fontfamily = "Helvetica",
         fontsize = 14, fontface = "bold")
plotGG(plot = skipped_locus_barplot, x = 4.75, y = 0.25, width = 4.25, height = 2.25)
dev.off()

