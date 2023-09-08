
# Sex-specific response to treatment -------------------------------------------
load("data/dds_sex_treatmentresponse.rda")

# Normalized counts
dds_sex_treatmentresponse_norm <- vst(dds_sex_treatmentresponse)

# Read in significant treated genes and order by log2FC and padj
treatmentresponse_sex_degenes <- read_csv("data/sexDE_treatmenteffect_pval05.csv") %>%
  arrange(log2FoldChange) %>%
  mutate(log2FC_dir = ifelse(log2FoldChange < 0, "-", "+")) %>%
  group_by(log2FC_dir) %>%
  arrange(padj, .by_group = TRUE)

# Subset norm counts for significant treated genes
sexnormCounts_treatmentresponse <- assay(dds_sex_treatmentresponse_norm[treatmentresponse_sex_degenes$gene_id,]) %>% 
  as.data.frame()
sexnormCounts_treatmentresponse <- sexnormCounts_treatmentresponse[match(treatmentresponse_sex_degenes$gene_id, 
                                                                         rownames(sexnormCounts_treatmentresponse)),]

# Reorder into M/F, CTL/FNF within each sex group
sexnormCounts_treatmentresponse_M <- sexnormCounts_treatmentresponse %>%
  dplyr::select(ends_with("_M")) %>%
  dplyr::select(contains(c("CTL", "FNF")))

sexnormCounts_treatmentresponse_F <- sexnormCounts_treatmentresponse %>%
  dplyr::select(ends_with("_F")) %>%
  dplyr::select(contains(c("CTL", "FNF")))


sexnormCounts_treatmentresponse <- bind_cols(sexnormCounts_treatmentresponse_M,
                                             sexnormCounts_treatmentresponse_F)


# Scale counts
treatmentresponse_sexmat_scaled <- t(apply(sexnormCounts_treatmentresponse, 1, scale))
colnames(treatmentresponse_sexmat_scaled) <- colnames(sexnormCounts_treatmentresponse)

# Age, Sex, and Race Clusters
annotations <- as.data.frame(colData(dds_sex_treatmentresponse)[,c("Condition", "Sex", "Age", "Race")]) 
# Put in same order as matrix
annotations <- annotations[match(colnames(treatmentresponse_sexmat_scaled), rownames(annotations)),]

ageColors <- ochre_pal(palette = "olsen_qual")(6)
raceColors <- sunset1(n = 5)
annotationObjects <- HeatmapAnnotation(
  df = annotations,
  col = list(Sex = c("M" = "#6F8CC7",
                     "F" = "#C06E8B"),
             Condition = c("CTL" = "#B8B8B8", 
                           "FNF" = "#4A4A4A"),
             Age = c("31_40" = ageColors[1],
                     "41_50" = ageColors[2],
                     "51_60" = ageColors[3],
                     "61_70" = ageColors[4],
                     "71_80" = ageColors[5],
                     "81_90" = ageColors[6]),
             Race = c("ASIAN" = raceColors[1], 
                      "BL" = raceColors[2],
                      "HISP" = raceColors[3], 
                      "C" = raceColors[4],
                      "ARAB" = raceColors[5],
                      "Unknown" = "grey")),
  annotation_name_gp = gpar(fontfamily = "Helvetica",
                            fontface = "bold",
                            fontsize = 8),
  which = "column"
)


h3 <- pheatmap(treatmentresponse_sexmat_scaled, 
               show_rownames = TRUE,
               top_annotation = annotationObjects,
               row_title = NULL,
               cluster_cols = FALSE,
               cluster_rows = TRUE,
               show_colnames = FALSE,
               color = colorRampPalette(c(
                 "#73B5F9",
                 "black",
                 "#F5D24D"))(7),
               breaks = seq(-3, 3),
               show_row_dend = FALSE,
               show_column_dend = FALSE,
               heatmap_legend_param = list(at = c(-3, 3),
                                           title = "Relative Expression",
                                           border = NA,
                                           title_gp = gpar(fontfamily = "Helvetica",
                                                           fontsize = 8,
                                                           fontface = "bold"),
                                           labels_gp = gpar(fontfamily = "Helvetica",
                                                            fontsize = 8),
                                           legend_height = unit(1.5, 'in'),
                                           grid_width = unit(0.125, "in")),
               column_title_gp = gpar(fontfamily = "Helvetica"))

draw(h3,  background = "transparent", show_annotation_legend = FALSE)



# Sex-specific response to treatment l2fc ---------------------------------
load("data/dds_sex_treatmentresponse_moregenefilter.rda")

treatmentresponse_sex_degenes <- read_csv("data/sexDE_treatmenteffect_moregenefilter_pval05.csv") %>%
  arrange(padj)

# Normalized counts
dds_sex_treatmentresponse_norm <- 
  as.data.frame(assay(vst(dds_sex_treatmentresponse_moregenefilter)[treatmentresponse_sex_degenes$gene_id,]))


# Get M and F log2FC between FNF and CTL
M_l2fc <- dds_sex_treatmentresponse_norm %>% 
  # Subset for Male counts
  dplyr::select(ends_with("_M")) %>%
  # Rename by pasting donor name and condition
  rename_with(~paste0(str_extract(., "AM\\d{4}"), "_", str_extract(., "CTL|FNF"))) %>%
  rownames_to_column(var = "gene_id") %>%
  # Long format, expanding donor and condition for each gene
  pivot_longer(cols = !gene_id,
               names_sep = "_",
               names_to = c("Donor", "Condition"),
               values_to = "count") %>%
  group_by(Donor, gene_id) %>%
  # Calculate log2fc between FNF and CTL for each donor of each gene
  summarize(log2FC = log2(count[Condition == "FNF"]/count[Condition == "CTL"])) %>%
  ungroup() %>% 
  # Wide format for plotting
  pivot_wider(names_from = Donor, 
              values_from = log2FC,
              names_glue = "{Donor}_M")

F_l2fc <- dds_sex_treatmentresponse_norm %>% 
  # Subset for Female counts
  dplyr::select(ends_with("_F")) %>%
  # Rename by pasting donor name and condition
  rename_with(~paste0(str_extract(., "AM\\d{4}"), "_", str_extract(., "CTL|FNF"))) %>%
  rownames_to_column(var = "gene_id") %>%
  # Long format, expanding donor and condition for each gene
  pivot_longer(cols = !gene_id,
               names_sep = "_",
               names_to = c("Donor", "Condition"),
               values_to = "count") %>%
  group_by(Donor, gene_id) %>%
  # Calculate log2fc between FNF and CTL for each donor of each gene
  summarize(log2FC = log2(count[Condition == "FNF"]/count[Condition == "CTL"])) %>%
  ungroup() %>% 
  # Wide format for plotting
  pivot_wider(names_from = Donor, 
              values_from = log2FC,
              names_glue = "{Donor}_F")


M_F_l2fc <- left_join(M_l2fc, F_l2fc, by = "gene_id") %>%
  column_to_rownames(var = "gene_id") %>%
  as.matrix()

# Order by padj
M_F_l2fc <- M_F_l2fc[match(treatmentresponse_sex_degenes$gene_id, rownames(M_F_l2fc)),]


annotations <- as.data.frame(colData(dds_sex_treatmentresponse_moregenefilter)[,c("Donor", "Condition", "Sex")]) %>%
  filter(Condition == "CTL") %>%
  mutate(id = paste0(Donor, "_", Sex)) %>%
  dplyr::select(-Condition,-Donor) %>%
  remove_rownames() %>%
  column_to_rownames(var = "id")

# Put in same order as matrix
annotations <- data.frame("Sex" = annotations[match(colnames(M_F_l2fc), 
                                                    rownames(annotations)),])


annotationObjects <- HeatmapAnnotation(
  df = annotations,
  col = list(Sex = c("M" = "#6F8CC7",
                     "F" = "#C06E8B")),
  annotation_name_gp = gpar(fontfamily = "Helvetica",
                            fontface = "bold",
                            fontsize = 8),
  which = "column"
)


h4 <- pheatmap(M_F_l2fc,
               show_rownames = FALSE,
               cluster_cols = FALSE,
               cluster_rows = FALSE,
               show_colnames = FALSE,
               show_row_dend = FALSE,
               show_column_dend = FALSE,
               color = colorRampPalette(c(
                 "#73B5F9",
                 "black",
                 "#F5D24D"))(3),
               breaks = seq(-1, 1),
               top_annotation = annotationObjects,
               heatmap_legend_param = list(at = c(-1, 1),
                                           title = "Sample FN-f log2FC",
                                           border = NA,
                                           title_gp = gpar(fontfamily = "Helvetica",
                                                           fontsize = 8,
                                                           fontface = "bold"),
                                           labels_gp = gpar(fontfamily = "Helvetica",
                                                            fontsize = 8),
                                           legend_height = unit(1.5, 'in'),
                                           grid_width = unit(0.125, "in")),
               column_title_gp = gpar(fontfamily = "Helvetica"))
draw(h4,  background = "transparent", show_annotation_legend = FALSE)