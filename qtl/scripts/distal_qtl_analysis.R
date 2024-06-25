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


# FUNCTIONS ---------------------------------------------------------------

# Function to connect a list of loops to genes 
# `expression` will get the corresponding differential expression data for those genes
connect_loops <- function(loops, loop_binsize, promoters_data, category,
                          expression = TRUE){
  
  # Adjust binsize of loops
  binsize_loops <- snapToBins(loops, binSize = loop_binsize)
  
  # Connect each loop to a gene
  loop_gene_overlaps <- findOverlaps(query = binsize_loops,
                                     subject = promoters_data)
  
  # Create version of loop dataset that have all the query hits
  loops_genes <- loops[queryHits(loop_gene_overlaps)]
  
  # Extract looped gene IDs and add as a column to loops_genes
  looped_genes <- unlist(promoters_data[subjectHits(loop_gene_overlaps)]$GENEID)
  loops_genes$gene_id <- looped_genes
  
  # Convert to tibble and get distinct loop_gene pairs
  loops_genes_df <- as_tibble(loops_genes) |> 
    mutate(loop_id = names(loops_genes)) |> 
    distinct()
  
  
  if (expression == TRUE){
    # Create separate df for plotting and join with RNA log2FC
    loop_gene_final <- tibble(loopCategory = category, 
                              loop_id = loops_genes_df$loop_id,
                              gene_id = loops_genes_df$gene_id) |> 
      # Get log2FC of genes from DE data
      left_join(read_csv("../DE/data/condition_de/de_genes_results.csv",
                         col_select = c("gene_id", "symbol", "log2FoldChange", "padj")),
                by = "gene_id") |> 
      # Filter ones we don't have DE info for
      filter(!is.na(log2FoldChange))
  } else {
    
    # Don't join with DE info 
    loop_gene_final <- tibble(loopCategory = category,
                              loop_id = loops_genes_df$loop_id,
                              gene_id = loops_genes_df$gene_id)
    
  }
  
  return(loop_gene_final)
  
}

# Function to categorize an independent eSNP-eGene signals as distal or not.
# This is defined with the signal range (variants in LD > 0.6) being at least
# 50 Kb away from either end of the gene.
determine_distal <- function(eSNP_eGene, txdb_genes_gr){
  eGene_id <- eSNP_eGene[["gene_id"]]
  signal_start <- as.numeric(eSNP_eGene[["min_signal_pos"]])
  signal_end <- as.numeric(eSNP_eGene[["max_signal_pos"]])
  
  # Grab eGene from txdb_genes_gr to get entire gene range
  eGene_range <- txdb_genes_gr |> 
    filter(gene_id == eGene_id)
  
  if (length(eGene_range) > 0){
    gene_upstream_50kb <- start(eGene_range) - 50000
    gene_downstream_50kb <- end(eGene_range) + 50000
    
    # Compare signal end to gene start
    if (signal_end <= gene_upstream_50kb | signal_start >= gene_downstream_50kb){
      return("distal")
    } else {
      return("not_distal")
    }
    
  } else {
    return("not_distal")
  }
  
}

# Function to check which eSNP-eGene connections are supported by loops
determine_loop_support <- function(eSNP_eGene, 
                                   promoters_mapped_genes, 
                                   loops,
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
  loop_interactions <- linkOverlaps(snapToBins(loops, binSize = loop_binsize), 
                                    subject1 = signal_range, 
                                    subject2 = eGene_promoters)
  
  if (nrow(loop_interactions) > 0){
    return("loop")
  } else {
    return("no_loop")
  }
  
}


# Function to check the loop support of individual signal SNPs to their eGene
determine_ld_loop_support <- function(eSNP_eGene, 
                                   promoters_mapped_genes, 
                                   loops,
                                   loop_binsize){
  eGene_id <- eSNP_eGene[["gene_id"]]
  variant_pos <- as.numeric(eSNP_eGene[["ld_pos"]])
  chrom <- eSNP_eGene[["variant_chr"]]
  
  # Make variant a GRanges object for overlapping
  eSNP_gr <- GRanges(seqnames = chrom,
                     ranges = IRanges(start = variant_pos,
                                      end = variant_pos))
  
  # Grab eGene_id promoter ranges
  eGene_promoters <- promoters_mapped_genes |> 
    filter(GENEID == eGene_id)
  
  # Get loop interactions that link snp and any of the eGene's promoters
  loop_interactions <- linkOverlaps(snapToBins(loops, binSize = loop_binsize), 
                                    subject1 = eSNP_gr, 
                                    subject2 = eGene_promoters)
  
  if (nrow(loop_interactions) > 0){
    return("loop")
  } else {
    return("no_loop")
  }
  
}

# Function to create a GInteractions object between a signal's
# lead variant and it's eGene's promoters 
make_lead_eGene_gi <- function(eSNP_eGene, promoters_mapped_genes){
  eGene_id <- eSNP_eGene[["gene_id"]]
  chrom <- eSNP_eGene[["gene_chr"]]
  lead_var_pos <- as.numeric(eSNP_eGene[["variant_start"]])
  variantID <- eSNP_eGene[["variantID"]]
  print(variantID)
  
  # Make GRanges for lead variant
  lead_var_gr <- GRanges(seqnames = chrom,
                         ranges = IRanges(start = lead_var_pos,
                                          end = lead_var_pos))
  
  
  # Grab eGene_id promoter ranges
  eGene_promoters <- promoters_mapped_genes |> 
    filter(GENEID == eGene_id)
  
  if (length(eGene_promoters) > 0){
    # Create GInteractions between the lead and these promoters
    lead_eGene_gr <- GInteractions(anchor1 = rep(lead_var_gr, length(eGene_promoters)), 
                                   anchor2 = eGene_promoters,
                                   variantID = variantID)
    # Rename anchor2.GENEID column and anchor2.TXID column
    lead_eGene_gr$gene_id <- lead_eGene_gr$anchor2.GENEID
    lead_eGene_gr$anchor2.GENEID <- NULL
    lead_eGene_gr$tx_id <- lead_eGene_gr$anchor2.TXID
    lead_eGene_gr$anchor2.TXID <- NULL
  } else {
    # If no gene info, create empty GInteractions
    lead_eGene_gr <- GInteractions(variantID = variantID, 
                                   gene_id = eGene_id, tx_id = eGene_id)
  }
  
  return(lead_eGene_gr)
  
}

# Function to create a GInteractions object between a signal's
# lead variant and other genes that aren't its eGene and are +- 1Mb away
make_lead_gene_gi <- function(eSNP_eGene, promoters_mapped_genes){
  eGene_id <- eSNP_eGene[["gene_id"]]
  chrom <- eSNP_eGene[["gene_chr"]]
  lead_var_pos <- as.numeric(eSNP_eGene[["variant_start"]])
  variantID <- eSNP_eGene[["variantID"]]
  
  # Make GRanges for lead variant
  lead_var_gr <- GRanges(seqnames = chrom,
                         ranges = IRanges(start = lead_var_pos,
                                          end = lead_var_pos))
  
  
  # Grab gene promoter ranges for every gene except the eGene
  gene_promoters <- promoters_mapped_genes |> 
    filter(seqnames == chrom & start > lead_var_pos - 1000000 & end < lead_var_pos + 1000000) |> 
    filter(GENEID != eGene_id)
  
  # Create GInteractions between the lead and these promoters
  lead_gene_gr <- GInteractions(anchor1 = rep(lead_var_gr, length(gene_promoters)), 
                                anchor2 = gene_promoters,
                                variantID = variantID)
  # Rename anchor2.GENEID column and anchor2.TXID column
  lead_gene_gr$gene_id <- lead_gene_gr$anchor2.GENEID
  lead_gene_gr$anchor2.GENEID <- NULL
  lead_gene_gr$tx_id <- lead_gene_gr$anchor2.TXID
  lead_gene_gr$anchor2.TXID <- NULL
  
  return(lead_gene_gr)
  
}

# Function to get any Hi-C loops that overlap an eQTL signal
get_signal_loops <- function(eQTL_signal, loops, loop_binsize){
  signal_start <- as.numeric(eQTL_signal[["min_signal_pos"]])
  signal_end <- as.numeric(eQTL_signal[["max_signal_pos"]])
  chrom <- eQTL_signal[["variant_chr"]]
  
  # Define signal range as GRanges object
  signal_range <- GRanges(seqnames = chrom,
                          ranges = IRanges(start = signal_start,
                                           end = signal_end))
  
  # Find loop overlaps with signal_range
  signal_loop_overlapHits <- findOverlaps(signal_range, snapToBins(loops, binSize = loop_binsize))
  signal_loop_overlaps <- loops[unique(subjectHits(signal_loop_overlapHits))]
  
  return(signal_loop_overlaps)
  
}

# RNA log2FC of genes at static loop anchors vs gained looped anchor --------

## Define Ensembl TxDb for reference
# Ensembl TxDb
ensembl_txdb <- makeTxDbFromEnsembl(organism = "Homo sapiens", release = 111)
# Change to UCSC style chromosomes
seqlevelsStyle(ensembl_txdb) <- "UCSC"
# Get promoters and filter for TXTYPE that contain any form of protein_coding
ensembl_promoters_coding <- promoters(ensembl_txdb, columns = c("TXID", "TXNAME", 
                                                                "GENEID", "TXTYPE")) |> 
  filter(grepl("*protein_coding*", TXTYPE))


## Load loops and add loop IDs based on rowname
# Gained
load("data/hic/CQTL_10kb_sig_gained_loops.rda")
names(CQTL_10kb_sig_gained_loops) <- paste0("gained_loop_", 
                                            seq(1:length(CQTL_10kb_sig_gained_loops)))

# Static
load("data/hic/CQTL_10kb_static_loops.rda")
names(CQTL_10kb_static_loops) <- paste0("static_loop_", 
                                        seq(1:length(CQTL_10kb_static_loops)))

## Call function to connect loops to gene promoters and get corresponding
## DE info

CQTL_10kb_gained_loop_expression <- connect_loops(loops = CQTL_10kb_sig_gained_loops,
                                                  loop_binsize = 30000,
                                                  promoters_data = ensembl_promoters_coding,
                                                  category = "gained")

CQTL_10kb_static_loop_expression <- connect_loops(loops = CQTL_10kb_static_loops,
                                                  loop_binsize = 30000,
                                                  promoters_data = ensembl_promoters_coding,
                                                  category = "static")

# Combine into one for plotting
CQTL_10kb_gained_static_loops_expression <- bind_rows(CQTL_10kb_gained_loop_expression,
                                                      CQTL_10kb_static_loop_expression)
CQTL_10kb_gained_static_loops_expression$loopCategory <- factor(CQTL_10kb_gained_static_loops_expression$loopCategory,
                                                                levels = c("static", "gained"))

write_csv(CQTL_10kb_gained_static_loops_expression, 
          "data/hic/CQTL_10kb_gained_static_loops_expression.csv")


# COMPILE LISTS OF GENES AT STATIC AND GAINED LOOP ANCHORS AND LAUNCH HOMER --------

## Define Ensembl TxDb for reference
# Ensembl TxDb
ensembl_txdb <- makeTxDbFromEnsembl(organism = "Homo sapiens", release = 110)
# Change to UCSC style chromosomes
seqlevelsStyle(ensembl_txdb) <- "UCSC"
# Get promoters and filter for TXTYPE that contain any form of protein_coding
ensembl_promoters_coding <- promoters(ensembl_txdb, columns = c("TXID", "TXNAME", 
                                                                "GENEID", "TXTYPE")) |> 
  filter(grepl("*protein_coding*", TXTYPE))

## Load loops and add loop IDs based on rowname
# Gained
load("data/hic/CQTL_10kb_sig_gained_loops.rda")
names(CQTL_10kb_sig_gained_loops) <- paste0("gained_loop_", seq(1:length(CQTL_10kb_sig_gained_loops)))

# Static
load("data/hic/CQTL_10kb_static_loops.rda")
names(CQTL_10kb_static_loops) <- paste0("static_loop_", seq(1:length(CQTL_10kb_static_loops)))

## Call function to connect loops to gene promoters and write to files for
## Homer gene list and background

genes_gained_loops <- connect_loops(loops = CQTL_10kb_sig_gained_loops,
                                    loop_binsize = 30000,
                                    promoters_data = ensembl_promoters_coding,
                                    category = "gained",
                                    expression = FALSE) |> 
  dplyr::select(gene_id) |> 
  write_csv(file = "data/hic/homer/genes_gained_loops.csv", col_names = FALSE)

genes_static_loops <- connect_loops(loops = CQTL_10kb_static_loops,
                                    loop_binsize = 30000,
                                    promoters_data = ensembl_promoters_coding,
                                    category = "static",
                                    expression = FALSE) |> 
  dplyr::select(gene_id) |> 
  write_csv(file = "data/hic/homer/genes_static_loops.csv", col_names = FALSE)

# Launch homer
system("scripts/gene_loopAnchors_homer.sh data/hic/homer/genes_gained_loops.csv data/hic/homer/genes_static_loops.csv data/hic/homer")

# DISTAL SIGNAL CLASSIFICATION AND LOOP SUPPORT --------------------------------

## Distal
# Read in pbs and fnf conditional signal data with calculated signal ranges
# generated from distal_Fig5_getSignalRanges.R
pbs_signals_signalRanges <- read_csv("data/eqtl/CTL_PEER_k20_genoPC_cond1Mb_topSignals_rsID_signalRanges.csv") |> 
  mutate(Condition = "PBS")
fnf_signals_signalRanges <- read_csv("data/eqtl/FNF_PEER_k22_genoPC_cond1Mb_topSignals_rsID_signalRanges.csv") |> 
  mutate(Condition = "FN-f")

# Get ENSEMBL Txdb
ensembl_txdb <- makeTxDbFromEnsembl(organism = "Homo sapiens", release = 111)
# Update seqlevels to UCSC style for compatibility with eqtl data
seqlevelsStyle(ensembl_txdb) <- "UCSC"
# Grab gene ranges
ensembl_genes <- genes(ensembl_txdb)

# Classify signals as distal if 50kb away from either end of eGene
pbs_signals_signalRanges$distal <- apply(pbs_signals_signalRanges, 1, 
                                         determine_distal, ensembl_genes)
fnf_signals_signalRanges$distal <- apply(fnf_signals_signalRanges, 1, 
                                         determine_distal, ensembl_genes)

# Write to files
write_csv(pbs_signals_signalRanges, 
          file = "data/eqtl/CTL_PEER_k20_genoPC_cond1Mb_topSignals_rsID_signalRanges_distal.csv")
write_csv(fnf_signals_signalRanges, 
          file = "data/eqtl/FNF_PEER_k22_genoPC_cond1Mb_topSignals_rsID_signalRanges_distal.csv")

# Join together
all_signals_signalRanges <- bind_rows(pbs_signals_signalRanges, 
                                      fnf_signals_signalRanges) |> 
  distinct() 


## Distal loop support
distal_all_signals_signalRanges <- all_signals_signalRanges  |> 
  filter(distal == "distal")

# Load static and sig gained loops and concat
load("data/hic/CQTL_10kb_static_loops.rda")
load("data/hic/CQTL_10kb_sig_gained_loops.rda")
CQTL_10kb_allLoops <- c(CQTL_10kb_static_loops, CQTL_10kb_sig_gained_loops)
# Add names to loops for easier access
names(CQTL_10kb_allLoops) <- c(paste0("static_loop_", 
                                      seq(1:length(CQTL_10kb_static_loops))),
                               paste0("gained_loop_", 
                                      seq(1:length(CQTL_10kb_sig_gained_loops))))

# Get ENSEMBL TxDb gene promoters
ensembl_promoter_genes <- promoters(ensembl_txdb, columns = c("GENEID", "TXID"))
# convert GENEID CharacterList just to character
ensembl_promoter_genes$GENEID <- as.character(ensembl_promoter_genes$GENEID)

# Add additional column for loop support with expanded 30kb anchors
distal_all_signals_signalRanges$looping <- apply(distal_all_signals_signalRanges, 
                                                 1, 
                                                 determine_loop_support, 
                                                 ensembl_promoter_genes,
                                                 CQTL_10kb_allLoops,
                                                 loop_binsize = 30000)
## Write to file
write_csv(distal_all_signals_signalRanges,
          file = "tables/eGenes_signals_distal_looping.csv")

# CONTACT FREQUENCY BETWEEN DISTAL LEAD SNPS AND THEIR EGENE PROMOTERS ---------  

#### TXDB PARSING
# Get ENSEMBL Txdb
ensembl_txdb <- makeTxDbFromEnsembl(organism = "Homo sapiens", release = 111)
# Only keep standard chromosomes for compatibility with Hi-C files
ensembl_txdb <- keepStandardChromosomes(ensembl_txdb)
# Update seqlevels to UCSC style for compatibility with eqtl data
seqlevelsStyle(ensembl_txdb) <- "UCSC"
# Get ENSEMBL TxDb gene promoters
ensembl_promoter_genes <- promoters(ensembl_txdb, columns = c("GENEID", "TXID"))
# convert GENEID CharacterList just to character
ensembl_promoter_genes$GENEID <- as.character(ensembl_promoter_genes$GENEID)

#### CREATE GINTERACTIONS OBJECTS BETWEEN DISTAL LEAD SNPS AND GENE PROMOTERS
#### FOR PIXEL EXTRACTION

# Create GInteractions between distal lead snps and eGene promoters
distal_all_signals_signalRanges <- read_csv("tables/eGenes_signals_distal_looping.csv")
lead_eGene_promoters <- apply(distal_all_signals_signalRanges, 1, make_lead_eGene_gi, ensembl_promoter_genes)
lead_eGene_promoters <- do.call(c, lead_eGene_promoters)

# Create GInteractions between distal lead snps and other genes within 2Mb window
# that aren't the eGene
lead_othergene_promoters <- apply(distal_all_signals_signalRanges, 1, make_lead_gene_gi, ensembl_promoter_genes)
lead_othergene_promoters <- do.call(c, lead_othergene_promoters)

# Make interactions 5kb
lead_eGene_promoters_5kb <- snapToBins(lead_eGene_promoters, binSize = 5000)
lead_othergene_promoters_5kb <- snapToBins(lead_othergene_promoters, binSize = 5000)

# Save objects
save(lead_eGene_promoters_5kb, file = "data/hic/lead_eGene_promoters_5kb.rda")
save(lead_othergene_promoters_5kb, file = "data/hic/lead_othergene_promoters_5kb.rda")

#### PULL PIXELS
lead_eGene_promoters_5kb_pixels <- pullHicPixels(lead_eGene_promoters_10kb, 
                                                 files = c("data/hic/CQTL_AM7682_AM7683_AM7697_AM7698_PBS_inter_30.hic",
                                                           "data/hic/CQTL_AM7682_AM7683_AM7697_AM7698_FNF_inter_30.hic"),
                                                 binSize = 5000, norm = "SCALE",
                                                 h5File = "data/hic/lead_eGene_promoter_5kb_pixels.h5")
save(lead_eGene_promoters_5kb_pixels, file = "data/hic/lead_eGene_promoter_5kb_pixels.rda")

lead_othergene_promoters_5kb_pixels <- pullHicPixels(lead_othergene_promoters_5kb, 
                                                     files = c("data/hic/CQTL_AM7682_AM7683_AM7697_AM7698_PBS_inter_30.hic",
                                                               "data/hic/CQTL_AM7682_AM7683_AM7697_AM7698_FNF_inter_30.hic"),
                                                     binSize = 5000, norm = "SCALE",
                                                     h5File = "data/hic/lead_othergene_promoter_5kb_pixels.h5")
save(lead_othergene_promoters_5kb_pixels, file = "data/hic/lead_othergene_promoter_5kb_pixels.rda")

#### COMBINE INTERACTIONS AND COUNT INFO

lead_eGene_promoters_5kb_pixels_counts <- interactions(lead_eGene_promoters_5kb_pixels) |> 
  as.data.frame() |> 
  bind_cols(counts(lead_eGene_promoters_5kb_pixels) |> as.data.frame()) |> 
  # Mark as being an eGene pixel
  mutate(gene_category = "eGene") |> 
  as_ginteractions() 


lead_othergene_promoters_5kb_pixels_counts <- interactions(lead_othergene_promoters_5kb_pixels) |> 
  as.data.frame() |> 
  bind_cols(counts(lead_othergene_promoters_5kb_pixels) |> as.data.frame()) |> 
  # Mark as being a pixel for another gene
  mutate(gene_category = "other") |> 
  as_ginteractions()

# Join both together
distal_leadVar_gene_pixel_counts <- c(lead_eGene_promoters_5kb_pixels_counts,
                                      lead_othergene_promoters_5kb_pixels_counts) 

# Calculate distances between variant range and gene promoter range 
distal_leadVar_gene_pixel_counts$ranges_distance <- calculateDistances(distal_leadVar_gene_pixel_counts)

# Calculate total contact frequency by summing counts across Hi-C files
distal_leadVar_gene_pixel_counts$contact_freq <- apply(mcols(distal_leadVar_gene_pixel_counts)[grep("*inter_30.hic",
                                                                                                    colnames(mcols(distal_leadVar_gene_pixel_counts)))], 1, sum)

# Condense to one gene
distal_leadVar_gene_pixel_counts <- distal_leadVar_gene_pixel_counts |> 
  as_tibble() |> 
  # Average contact frequency and range distance across variant/gene
  group_by(gene_id, variantID) |> 
  mutate(mean_contact_freq = mean(contact_freq),
         mean_range_distance = mean(ranges_distance)) |> 
  ungroup() |> 
  # Keep first instance of gene
  distinct(variantID, gene_id, .keep_all = TRUE) |> 
  as.data.frame() |> 
  as_ginteractions()

#### NULLRANGES FOR DISTANCE-MATCHING
# Generate pixel matches based on the ranges_distance between variant and gene
# Focal is eGene pixels and pool is other gene pixels
distal_leadVar_pixel_matches <- matchRanges(focal = distal_leadVar_gene_pixel_counts[distal_leadVar_gene_pixel_counts$gene_category == "eGene"],
                                            pool = distal_leadVar_gene_pixel_counts[distal_leadVar_gene_pixel_counts$gene_category == "other"],
                                            covar = ~mean_range_distance, replace = FALSE, method = "rejection")

save(distal_leadVar_pixel_matches, 
     file = "data/hic/distal_leadVar_pixel_matches.rda")

# CONTACT FREQUENCY OF DISTAL PBS-SPECIFIC, SHARED, AND FNF-SPECIFIC EQTLS -----

distal_all_signals_signalRanges <- read_csv("tables/eGenes_signals_distal_looping.csv")

#### TXDB PARSING
# Get ENSEMBL Txdb
ensembl_txdb <- makeTxDbFromEnsembl(organism = "Homo sapiens", release = 111)
# Only keep standard chromosomes for compatibility with Hi-C files
ensembl_txdb <- keepStandardChromosomes(ensembl_txdb)
# Update seqlevels to UCSC style for compatibility with eqtl data
seqlevelsStyle(ensembl_txdb) <- "UCSC"
# Get ENSEMBL TxDb gene promoters
ensembl_promoter_genes <- promoters(ensembl_txdb, columns = c("GENEID", "TXID"))
# convert GENEID CharacterList just to character
ensembl_promoter_genes$GENEID <- as.character(ensembl_promoter_genes$GENEID)

#### Distal high-confidence PBS and FN-f-specific eGenes
pbs_distal_highconf_regenes <- read_csv("data/reqtl/CTL_sig01_beta_donor_reQTLs_PEER_k20_genoPC.csv") |> 
  filter(gene_id %in% distal_all_signals_signalRanges$gene_id)
fnf_distal_highconf_regenes <- read_csv("data/reqtl/FNF_sig01_beta_donor_reQTLs_PEER_k22_genoPC.csv") |> 
  filter(gene_id %in% distal_all_signals_signalRanges$gene_id)

# Distal shared eGenes
pbs_distal_egenes <- read_csv("data/eqtl/CTL_PEER_k20_genoPC_perm1Mb_sig_rsID.csv") |> 
  filter(gene_id %in% distal_all_signals_signalRanges$gene_id) |> 
  pull(gene_id)
fnf_distal_egenes <- read_csv("data/eqtl/FNF_PEER_k22_genoPC_perm1Mb_sig_rsID.csv") |> 
  filter(gene_id %in% distal_all_signals_signalRanges$gene_id) |> 
  pull(gene_id)

shared_distal_egenes <- intersect(pbs_distal_egenes, fnf_distal_egenes)

## Get lead eQTLs and LD > 0.8 for each category
# PBS
highconf_PBS_distal_signals_LD <- fread("data/eqtl/CTL_PEER_k20_genoPC_cond1Mb_topSignals_rsID_LD.csv", 
                                        data.table = FALSE) |> 
  filter(gene_id %in% pbs_distal_highconf_regenes$gene_id) |> 
  filter(R2 > 0.8) |> 
  separate_wider_delim(cols = "ld_variantID", delim = ":", names = c(NA, "ld_pos", NA, NA), cols_remove = FALSE) |> 
  dplyr::select(gene_id, gene_chr, ld_pos, ld_variantID) |> 
  dplyr::rename(variant_start = ld_pos,
                variantID = ld_variantID)

# FNF
highconf_FNF_distal_signals_LD <- fread("data/eqtl/FNF_PEER_k22_genoPC_cond1Mb_topSignals_rsID_LD.csv",
                                        data.table = FALSE) |> 
  filter(gene_id %in% fnf_distal_highconf_regenes$gene_id) |> 
  filter(R2 > 0.8) |> 
  separate_wider_delim(cols = "ld_variantID", delim = ":", names = c(NA, "ld_pos", NA, NA), cols_remove = FALSE) |> 
  dplyr::select(gene_id, gene_chr, ld_pos, ld_variantID) |> 
  dplyr::rename(variant_start = ld_pos,
                variantID = ld_variantID)

# Shared from either PBS or FNF
shared_distal_signals_LD <- bind_rows(fread("data/eqtl/CTL_PEER_k20_genoPC_cond1Mb_topSignals_rsID_LD.csv", 
                                            data.table = FALSE) |> 
                                        filter(gene_id %in% shared_distal_egenes) |> 
                                        filter(R2 > 0.8) |> 
                                        separate_wider_delim(cols = "ld_variantID", 
                                                             delim = ":",
                                                             names = c(NA, "ld_pos", NA, NA), cols_remove = FALSE) |> 
                                        dplyr::select(gene_id, ld_variantID, gene_chr, ld_pos),
                                      fread("data/eqtl/FNF_PEER_k22_genoPC_cond1Mb_topSignals_rsID_LD.csv", 
                                            data.table = FALSE) |> 
                                        filter(gene_id %in% shared_distal_egenes) |>
                                        filter(R2 > 0.8) |> 
                                        separate_wider_delim(cols = "ld_variantID", 
                                                             delim = ":",
                                                             names = c(NA, "ld_pos", NA, NA), cols_remove = FALSE) |>
                                        dplyr::select(gene_id, ld_variantID, gene_chr, ld_pos)) |> 
  distinct() |> 
  dplyr::rename(variant_start = ld_pos,
                variantID = ld_variantID)

## Create GInteractions between snps and eGene promoters and make 5kb for pulling pixels
highconf_PBS_distal_signals_LD_eGenes <- apply(highconf_PBS_distal_signals_LD, 1, make_lead_eGene_gi, ensembl_promoter_genes)
highconf_PBS_distal_signals_LD_eGenes <- do.call(c, highconf_PBS_distal_signals_LD_eGenes)
highconf_PBS_distal_signals_LD_eGenes_5kb <- snapToBins(highconf_PBS_distal_signals_LD_eGenes,
                                                        binSize = 5000)
save(highconf_PBS_distal_signals_LD_eGenes_5kb, file = "data/hic/highconf_PBS_distal_signals_LD_eGenes_5kb.rda")

highconf_FNF_distal_signals_LD_eGenes <- apply(highconf_FNF_distal_signals_LD, 1, make_lead_eGene_gi, ensembl_promoter_genes)
highconf_FNF_distal_signals_LD_eGenes <- do.call(c, highconf_FNF_distal_signals_LD_eGenes)
highconf_FNF_distal_signals_LD_eGenes_5kb <- snapToBins(highconf_FNF_distal_signals_LD_eGenes, 
                                                        binSize = 5000)
save(highconf_FNF_distal_signals_LD_eGenes_5kb, file = "data/hic/highconf_FNF_distal_signals_LD_eGenes_5kb.rda")

shared_distal_signals_LD_eGenes <- apply(shared_distal_signals_LD, 1, make_lead_eGene_gi, ensembl_promoter_genes)
shared_distal_signals_LD_eGenes <- do.call(c, shared_distal_signals_LD_eGenes)
shared_distal_signals_LD_eGenes_5kb <- snapToBins(shared_distal_signals_LD_eGenes,
                                                  binSize = 5000)
shared_distal_signals_LD_eGenes_5kb[width(shared_distal_signals_LD_eGenes_5kb)$first == 1] <- binPairs(shared_distal_signals_LD_eGenes_5kb[width(shared_distal_signals_LD_eGenes_5kb)$first == 1], binSize = 5000)
save(shared_distal_signals_LD_eGenes_5kb, file = "data/hic/shared_distal_signals_LD_eGenes_5kb.rda")

## Pull pixels for all these interactions at 5kb
highconf_PBS_distal_signals_LD_eGenes_pixels_5kb <- pullHicPixels(highconf_PBS_distal_signals_LD_eGenes_5kb, 
                                                                  files = c("data/hic/CQTL_AM7682_AM7683_AM7697_AM7698_PBS_inter_30.hic",
                                                                            "data/hic/CQTL_AM7682_AM7683_AM7697_AM7698_FNF_inter_30.hic"),
                                                                  binSize = 5000, norm = "SCALE",
                                                                  h5File = "data/hic/highconf_PBS_distal_signals_LD_eGenes_pixels_5kb.h5")
save(highconf_PBS_distal_signals_LD_eGenes_pixels_5kb, file = "data/hic/highconf_PBS_distal_signals_LD_eGenes_pixels_5kb.rda")

highconf_FNF_distal_signals_LD_eGenes_pixels_5kb <- pullHicPixels(highconf_FNF_distal_signals_LD_eGenes_5kb, 
                                                                  files = c("data/hic/CQTL_AM7682_AM7683_AM7697_AM7698_PBS_inter_30.hic",
                                                                            "data/hic/CQTL_AM7682_AM7683_AM7697_AM7698_FNF_inter_30.hic"),
                                                                  binSize = 5000, norm = "SCALE",
                                                                  h5File = "data/hic/highconf_FNF_distal_signals_LD_eGenes_pixels_5kb.h5")
save(highconf_FNF_distal_signals_LD_eGenes_pixels_5kb, file = "data/hic/highconf_FNF_distal_signals_LD_eGenes_pixels_5kb.rda")

shared_distal_signals_LD_eGenes_pixels_5kb <- pullHicPixels(shared_distal_signals_LD_eGenes_5kb, 
                                                            files = c("data/hic/CQTL_AM7682_AM7683_AM7697_AM7698_PBS_inter_30.hic",
                                                                      "data/hic/CQTL_AM7682_AM7683_AM7697_AM7698_FNF_inter_30.hic"),
                                                            binSize = 5000, norm = "SCALE",
                                                            h5File = "data/hic/shared_distal_signals_LD_eGenes_pixels_5kb.h5")
save(shared_distal_signals_LD_eGenes_pixels_5kb, file = "data/hic/shared_distal_signals_LD_eGenes_pixels_5kb.rda")

## Combine interactions and count info for all QTL categories
highconf_PBS_distal_signals_LD_eGenes_pixels_counts <- interactions(highconf_PBS_distal_signals_LD_eGenes_pixels_5kb) |> 
  as.data.frame() |> 
  bind_cols(counts(highconf_PBS_distal_signals_LD_eGenes_pixels_5kb) |> as.data.frame()) |> 
  # Mark as being a high conf PBS-specific pixel
  mutate(eqtl_category = "PBS-specific")

highconf_FNF_distal_signals_LD_eGenes_pixels_counts <- interactions(highconf_FNF_distal_signals_LD_eGenes_pixels_5kb) |> 
  as.data.frame() |> 
  bind_cols(counts(highconf_FNF_distal_signals_LD_eGenes_pixels_5kb) |> as.data.frame()) |> 
  # Mark as being a high conf FNF-specific pixel
  mutate(eqtl_category = "FNF-specific")

shared_distal_signals_LD_eGenes_pixels_counts <- interactions(shared_distal_signals_LD_eGenes_pixels_5kb) |> 
  as.data.frame() |> 
  bind_cols(counts(shared_distal_signals_LD_eGenes_pixels_5kb) |> as.data.frame()) |> 
  # Mark as being a shared eqtl pixel
  mutate(eqtl_category = "shared")

## Join all datasets together
eqtl_categories_distal_egene_pixels_counts <- bind_rows(highconf_PBS_distal_signals_LD_eGenes_pixels_counts,
                                                        highconf_FNF_distal_signals_LD_eGenes_pixels_counts,
                                                        shared_distal_signals_LD_eGenes_pixels_counts) |> 
  mutate(across(ends_with("inter_30.hic"), as.numeric)) |> 
  # Calculate log2FC contact frequency between FNF and PBS Hi-C files
  mutate(log2FC_contact_freq = 
           log2(`CQTL_AM7682_AM7683_AM7697_AM7698_FNF_inter_30.hic`/`CQTL_AM7682_AM7683_AM7697_AM7698_PBS_inter_30.hic`)) |>
  # Filter out any Inf log2FC contact freq
  filter(is.finite(log2FC_contact_freq)) 

## Condense data down from multiple eQTL-gene promoter connections to
## one eQTL-gene connection
eqtl_categories_distal_egene_pixels_counts <- eqtl_categories_distal_egene_pixels_counts |> 
  # Average log2FC contact frequency across variant/gene
  group_by(gene_id, variantID) |> 
  mutate(mean_log2FC_contact_freq = mean(log2FC_contact_freq)) |> 
  ungroup() |> 
  # Keep first instance of gene
  distinct(variantID, gene_id, .keep_all = TRUE) |> 
  # Make group names based on eqtl_category
  mutate(eqtl_category_label = ifelse(eqtl_category == "PBS-specific", "**PBS-specific**<br> distal eQTLs<br> and eGenes",
                                      ifelse(eqtl_category == "FNF-specific", "**FN-f-specific**<br> distal eQTLs<br> and eGenes", 
                                             "**shared**<br> distal eQTLs<br> and eGenes"))) |> 
  # Change to factor correct plotting order
  mutate(eqtl_category_label = factor(eqtl_category_label, 
                                      levels = c("**PBS-specific**<br> distal eQTLs<br> and eGenes", 
                                                 "**shared**<br> distal eQTLs<br> and eGenes",
                                                 "**FN-f-specific**<br> distal eQTLs<br> and eGenes")))

write_csv(eqtl_categories_distal_egene_pixels_counts, file =
            "data/hic/eqtl_categories_distal_egene_pixels_counts.csv")

# DISTAL LOOP CONTACT FREQUENCY -------------------------------------------

# Get ENSEMBL Txdb
ensembl_txdb <- makeTxDbFromEnsembl(organism = "Homo sapiens", release = 111)
# Update seqlevels to UCSC style for compatibility with eqtl data
seqlevelsStyle(ensembl_txdb) <- "UCSC"
# Get ENSEMBL TxDb gene promoters
ensembl_promoter_genes <- promoters(ensembl_txdb, columns = c("GENEID", "TXID"))
# convert GENEID CharacterList just to character
ensembl_promoter_genes$GENEID <- as.character(ensembl_promoter_genes$GENEID)

# Get loops for each distal eSNP with loop anchors expanded to 30Kb
distal_egene_10kb_allLoops <- apply(distal_all_signals_signalRanges, 1, get_signal_loops, 
                                    loops = CQTL_10kb_allLoops,
                                    loop_binsize = 30000)
# Combine and get unique set
distal_egene_10kb_allLoops <- do.call(c, distal_egene_10kb_allLoops) |> 
  unique() 



# Get the genes these loops are connected to
distal_egene_10kb_allLoops_genes <- connect_loops(distal_egene_10kb_allLoops, 
                                                  loop_binsize = 30000,
                                                  promoters_data = ensembl_promoter_genes,
                                                  category = "distal_egenes",
                                                  expression = FALSE) |> 
  dplyr::select(-loopCategory)

# Subset these loops for ones also connected to a gene
distal_egene_10kb_loops_with_genes <- distal_egene_10kb_allLoops[names(distal_egene_10kb_allLoops) %in%
                                                                   distal_egene_10kb_allLoops_genes$loop_id] 

# Sum counts across samples and set as contact frequencey
distal_egene_10kb_loops_with_genes$contact_freq <- apply(mcols(distal_egene_10kb_loops_with_genes)[grep("*inter_30.hic",
                                                                                                        colnames(mcols(distal_egene_10kb_loops_with_genes)))], 1, sum)

# Mark which ones are signal_eGene loops
signal_eGene_loops <- apply(distal_all_signals_signalRanges, 1, get_signal_eGene_loops,
                            promoters_mapped_genes = ensembl_promoter_genes,
                            loops = CQTL_10kb_allLoops,
                            loop_binsize = 30000)
# Collapse into single GInteractions object
signal_eGene_loops <- do.call(c, do.call(c, lapply(signal_eGene_loops, 
                                                   function(x) x[!is.na(x)])))

distal_egene_10kb_loops_with_genes$signal_eGene_loop <- names(distal_egene_10kb_loops_with_genes) %in% names(signal_eGene_loops)

# Calculate loop distance
distal_egene_10kb_loops_with_genes$loop_distance <- start2(distal_egene_10kb_loops_with_genes) - start1(distal_egene_10kb_loops_with_genes)

# Generate matches based on loop distance  
distal_loop_matches <- matchRanges(focal = distal_egene_10kb_loops_with_genes[distal_egene_10kb_loops_with_genes$signal_eGene_loop],
                                   pool = distal_egene_10kb_loops_with_genes[!distal_egene_10kb_loops_with_genes$signal_eGene_loop],
                                   covar = ~loop_distance)  

# Grab data for plotting boxplots
distal_looping_contactfreq <- c(focal(distal_loop_matches), distal_loop_matches) |> 
  # Make tibble for tidy manipulation
  as_tibble() |> 
  # Add back loop IDs
  mutate(loop_id = c(names(focal(distal_loop_matches)), names(distal_loop_matches))) |> 
  # Make loop_id the first column
  relocate(loop_id) |> 
  # Check for distinct loops-gene connections
  distinct() |> 
  # Grab relevant columns (loop_id, sample counts(contact freq), and whether it was a signal_eGene_loop)
  dplyr::select(loop_id, contact_freq, signal_eGene_loop) |> 
  mutate(log2FC_contact_freq = log2(contact_freq)) |> 
  # Change group names of signal_eGene_loop 
  mutate(signal_eGene_loop = ifelse(signal_eGene_loop, "loop between<br> signal<br> and eGene", 
                                    "loop between<br> signal and<br> distance-matched<br> gene")) |> 
  # Change to factor correct plotting order
  mutate(signal_eGene_loop = factor(signal_eGene_loop, 
                                    levels = c("loop between<br> signal<br> and eGene", 
                                               "loop between<br> signal and<br> distance-matched<br> gene")))

# Calculate p-value
wilcox_contact_freq_test <- wilcox.test(x = distal_looping_contactfreq |> 
                                          filter(signal_eGene_loop == "loop between<br> signal<br> and eGene") |> 
                                          pull(log2FC_contact_freq),
                                        y = distal_looping_contactfreq |>  
                                          filter(signal_eGene_loop == "loop between<br> signal and<br> distance-matched<br> gene") |> 
                                          pull(log2FC_contact_freq),
                                        alternative = "greater")


# How many loops connect any signal SNPs with LD > 0.6 to eGene promoter? --------

# PBS eQTLs and LD > 0.8
pbs_eqtls_LD <- fread("data/eqtl/CTL_PEER_k20_genoPC_cond1Mb_topSignals_rsID_LD.csv",
                      data.table = FALSE) |> 
  filter(R2 > 0.6) |> 
  separate_wider_delim(cols = "ld_variantID", delim = ":", names = c(NA, "ld_pos", NA, NA), cols_remove = FALSE)

# FNF eQTLs and LD > 0.8
fnf_eqtls_LD <- fread("data/eqtl/FNF_PEER_k22_genoPC_cond1Mb_topSignals_rsID_LD.csv",
                      data.table = FALSE) |> 
  filter(R2 > 0.6) |> 
  separate_wider_delim(cols = "ld_variantID", delim = ":", names = c(NA, "ld_pos", NA, NA), cols_remove = FALSE)

# Bind PBS and FN-f
PBS_FNF_eqtls_LD <- bind_rows(pbs_eqtls_LD, fnf_eqtls_LD) |> 
  distinct(gene_id, ld_variantID, .keep_all = TRUE)

# Load static and sig gained loops and concat
load("data/hic/CQTL_10kb_static_loops.rda")
load("data/hic/CQTL_10kb_sig_gained_loops.rda")
CQTL_10kb_allLoops <- c(CQTL_10kb_static_loops, CQTL_10kb_sig_gained_loops)
# Add names to loops for easier access
names(CQTL_10kb_allLoops) <- c(paste0("static_loop_", 
                                      seq(1:length(CQTL_10kb_static_loops))),
                               paste0("gained_loop_", 
                                      seq(1:length(CQTL_10kb_sig_gained_loops))))

# Get ENSEMBL TxDb gene promoters
ensembl_txdb <- makeTxDbFromEnsembl(organism = "Homo sapiens", release = 111)
# Update seqlevels to UCSC style for compatibility with eqtl data
seqlevelsStyle(ensembl_txdb) <- "UCSC"
ensembl_promoter_genes <- promoters(ensembl_txdb, columns = c("GENEID", "TXID"))
# convert GENEID CharacterList just to character
ensembl_promoter_genes$GENEID <- as.character(ensembl_promoter_genes$GENEID)

# Add additional column for loop support
PBS_FNF_eqtls_LD$looping <- apply(PBS_FNF_eqtls_LD, 
                                  1, 
                                  determine_ld_loop_support, 
                                  ensembl_promoter_genes,
                                  CQTL_10kb_allLoops,
                                  loop_binsize = 10000)


