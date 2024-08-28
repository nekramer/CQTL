library(tidyverse)
library(mariner)
library(nullranges)
library(InteractionSet)
library(GenomicInteractions)

load("data/hic/lead_ld_eGene_promoter_5kb_pixels.rda")
load("data/hic/lead_ld_othergene_promoter_5kb_pixels.rda")


#### COMBINE INTERACTIONS AND COUNT INFO

lead_ld_eGene_promoters_5kb_pixels_counts <- interactions(lead_ld_eGene_promoters_5kb_pixels) |> 
  as.data.frame() |> 
  bind_cols(counts(lead_ld_eGene_promoters_5kb_pixels) |> as.data.frame()) |> 
  # Mark as being an eGene pixel
  mutate(gene_category = "eGene") |> 
  as_ginteractions() 
print('here1')

lead_ld_othergene_promoters_5kb_pixels_counts <- interactions(lead_ld_othergene_promoters_5kb_pixels) |> 
  as.data.frame() |> 
  bind_cols(counts(lead_ld_othergene_promoters_5kb_pixels) |> as.data.frame()) |> 
  # Mark as being a pixel for another gene
  mutate(gene_category = "other") |> 
  as_ginteractions()

print('here2')


# Join both together
distal_leadVar_ld_gene_pixel_counts <- c(lead_ld_eGene_promoters_5kb_pixels_counts,
                                         lead_ld_othergene_promoters_5kb_pixels_counts) 

print("here3")
# Calculate distances between variant range and gene promoter range 
distal_leadVar_ld_gene_pixel_counts$ranges_distance <- calculateDistances(distal_leadVar_ld_gene_pixel_counts)
print("here4")
# Calculate total contact frequency by summing counts across Hi-C files
distal_leadVar_ld_gene_pixel_counts$contact_freq <- apply(mcols(distal_leadVar_ld_gene_pixel_counts)[grep("*inter_30.hic",
                                                                                                          colnames(mcols(distal_leadVar_ld_gene_pixel_counts)))], 1, sum)
print("here5")

# Condense to one gene
distal_leadVar_ld_gene_pixel_counts <- distal_leadVar_ld_gene_pixel_counts |> 
  as_tibble() |> 
  # Average contact frequency and range distance across variant/gene
  group_by(gene_id, ld_variantID) |> 
  mutate(mean_contact_freq = mean(contact_freq),
         mean_range_distance = mean(ranges_distance)) |> 
  ungroup() |> 
  # Keep first instance of gene
  distinct(ld_variantID, gene_id, .keep_all = TRUE) |> 
  as.data.frame() |> 
  as_ginteractions()

print("here6")
distal_leadVar_ld_pixel_matches <- matchRanges(focal = distal_leadVar_ld_gene_pixel_counts[distal_leadVar_ld_gene_pixel_counts$gene_category == "eGene"],
                                               pool = distal_leadVar_ld_gene_pixel_counts[distal_leadVar_ld_gene_pixel_counts$gene_category == "other"],
                                               covar = ~mean_range_distance, replace = FALSE, method = "rejection")

save(distal_leadVar_ld_pixel_matches, 
     file = "data/hic/distal_leadVar_ld_pixel_matches.rda")