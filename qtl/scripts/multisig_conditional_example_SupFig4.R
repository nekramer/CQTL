library(tidyverse)
library(data.table)
library(plotgardener)


## MMP16
mmp16_cond <- read_csv("data/eqtl/FNF_PEER_k22_genoPC_cond1Mb_topSignals_rsID.csv") |> 
  filter(gene_symbol == "MMP16")

mmp16_chrom <- mmp16_cond |> pull(gene_chr) |> unique()
mmp16_gene_id <- mmp16_cond |> pull(gene_id) |> unique()
mmp16_sig1_topVar <- mmp16_cond |> filter(signal == 0) |> pull(variantID)
mmp16_sig2_topVar <- mmp16_cond |> filter(signal == 1) |> pull(variantID)
mmp16_sig1_topVar_rsID <- mmp16_cond |> filter(signal == 0) |> pull(rsID)
mmp16_sig2_topVar_rsID <- mmp16_cond |> filter(signal == 1) |> pull(rsID)
mmp16_gene_tss <- mmp16_cond |> pull(gene_start) |> unique()

## LD for both top MMP16 top variants
mmp16_cond_LD <- fread("data/eqtl/FNF_PEER_k22_genoPC_cond1Mb_topSignals_rsID_LD.csv",
                       data.table = FALSE) |> 
  filter(gene_symbol == "MMP16") |> 
  dplyr::select(signal, ld_variantID, R2) |> 
  dplyr::rename(variantID = ld_variantID)

######## Nominal data before conditional analysis
mmp16_fnf_nom_noCond <- fread(paste0("data/eqtl/qtl_nom/FNF_PEER_k22_genoPC_nom1Mb_", 
                                     mmp16_chrom, ".csv"), data.table = FALSE) |> 
  filter(gene_id == mmp16_gene_id) |> 
  # Join with LD for first signal
  left_join(mmp16_cond_LD |> filter(signal == 0), by = "variantID") |> 
  # Add rsIDs for top variants
  mutate(rsID = ifelse(variantID == mmp16_sig1_topVar, mmp16_sig1_topVar_rsID,
                       ifelse(variantID == mmp16_sig2_topVar, mmp16_sig2_topVar_rsID, NA))) |> 
  # Rename columns for pg
  dplyr::rename(snp = rsID,
                chrom = variant_chr,
                pos = variant_start,
                p = nom_pval) |> 
  # LD groups
  mutate(LDgrp = cut(R2, c(0, 0.2, 0.4, 0.6, 0.8, 1)))
mmp16_fnf_nom_noCond$LDgrp <- addNA(mmp16_fnf_nom_noCond$LDgrp)

####### Nominal data after conditional analysis
mmp16_fnf_nom_Cond <- fread(paste0("data/eqtl/qtl_cond/FNF_PEER_k22_genoPC_allSignals_nom1Mb_MAFs_", 
                                   mmp16_chrom, ".csv"), data.table = FALSE) |> 
  filter(gene_id == mmp16_gene_id) |> 
  # Join with LD for both signals 
  left_join(mmp16_cond_LD, by = c("signal", "variantID")) |> 
  # Add rsIDs for top variants
  mutate(rsID = ifelse(variantID == mmp16_sig1_topVar, mmp16_sig1_topVar_rsID,
                       ifelse(variantID == mmp16_sig2_topVar, mmp16_sig2_topVar_rsID, NA))) |> 
  # Rename columns for pg
  dplyr::rename(snp = rsID,
                chrom = variant_chr,
                pos = variant_start,
                p = nom_pval) |> 
  # LD groups
  mutate(LDgrp = cut(R2, c(0, 0.2, 0.4, 0.6, 0.8, 1)))
mmp16_fnf_nom_Cond$LDgrp <- addNA(mmp16_fnf_nom_Cond$LDgrp)

## eGene nominal threshold
mmp16_nom_threshold <- mmp16_fnf_nom_Cond |> pull(pval_nominal_threshold) |> unique()

####### Plotting all signals
mmp16_signal_region <- pgParams(chrom = mmp16_chrom, 
                                chromstart = mmp16_gene_tss - 500000,
                                chromend = mmp16_gene_tss + 500000,
                                x = 0.5, width = 4)

mmp16_signal_ylim <- ceiling(max(log10(c(mmp16_fnf_nom_Cond$p, mmp16_fnf_nom_noCond$p))*-1)) + 2

pdf("plots/reQTLs_Fig3/SupFig4.pdf", width = 5, height = 5.6)
pageCreate(width = 5, height = 5.6, showGuides = FALSE)

# Before conditional
plotText(label = "FN-f eQTL signal before conditional analysis",
         x = 2.5, y = 0.5, fontfamily = "Helvetica",
         fontsize = 10)
man_noCond <- plotManhattan(data = mmp16_fnf_nom_noCond,
                            params = mmp16_signal_region,
                            fill = colorby("LDgrp",
                              palette = colorRampPalette(c("#262C74",
                                                           "#98CDED",
                                                           "#499A53",
                                                           "#EEA741",
                                                           "#DD3931",
                                                           "grey"))),
                            range = c(0, mmp16_signal_ylim),
                            y = 0.5, height = 1.25,
                            snpHighlights = data.frame(snp = c(mmp16_sig1_topVar_rsID, mmp16_sig2_topVar_rsID),
                                                       pch = c(23, 23), 
                                                       cex = c(0.5, 0.5),
                                                       col = c("#DD3931", "#DD3931")))

annoYaxis(plot = man_noCond, at = seq(0, mmp16_signal_ylim, 2), 
          axisLine = TRUE, fontsize = 8)
plotText(
  label = "-log10(p-value)", x = 0.2, y = 1.2, rot = 90,
  fontsize = 6, just = "center",
  default.units = "inches"
)

plotText(label = mmp16_sig1_topVar_rsID, x = 1.8, y = 0.9,
         fontfamily = "Helvetica", fontsize = 7)
plotText(label = mmp16_sig2_topVar_rsID, x = 3.1, y = 1.1,
         fontfamily = "Helvetica", fontsize = 7)

plotLegend(legend = c("0.8 - 1.0",
                      "0.6 - 0.8",
                      "0.4 - 0.6",
                      "0.2 - 0.4",
                      "0.0 - 0.2"),
           fill = c("#DD3931", "#EEA741", "#499A53","#98CDED","#262C74"),
           x = 3.6, y = 0.85, width = 0.1, height = 0.4, border = FALSE, 
           fontsize = 6)
plotText(paste0("r2 relative to ", mmp16_sig1_topVar_rsID), x = 4.2, y = 0.75, fontsize = 6, fontfamily = "Helvetica")

# After conditional, signal 1
plotText(label = paste0("FN-f independent eQTL signal 1 conditioned on ", mmp16_sig2_topVar_rsID),
         x = 2.5, y = 1.9, fontfamily = "Helvetica",
         fontsize = 10)
man_condSig1 <- plotManhattan(data = mmp16_fnf_nom_Cond |> filter(signal == 0),
                            params = mmp16_signal_region,
                            fill = colorby("LDgrp",
                                           palette = colorRampPalette(c("#262C74",
                                                                        "#98CDED",
                                                                        "#499A53",
                                                                        "#EEA741",
                                                                        "#DD3931",
                                                                        "grey"))),
                            range = c(0, mmp16_signal_ylim),
                            y = 2, height = 1.25,
                            snpHighlights = data.frame(snp = c(mmp16_sig1_topVar_rsID, mmp16_sig2_topVar_rsID),
                                                       pch = c(23, 23), 
                                                       cex = c(0.5, 0.5),
                                                       col = c("#DD3931", "#DD3931")))

annoYaxis(plot = man_condSig1, at = seq(0, mmp16_signal_ylim, 2), 
          axisLine = TRUE, fontsize = 8)
plotText(
  label = "-log10(p-value)", x = 0.2, y = 2.7, rot = 90,
  fontsize = 6, just = "center",
  default.units = "inches"
)
plotText(label = mmp16_sig1_topVar_rsID, x = 1.8, y = 2.2,
         fontfamily = "Helvetica", fontsize = 7)
plotText(label = mmp16_sig2_topVar_rsID, x = 3.1, y = 3.1,
         fontfamily = "Helvetica", fontsize = 7)


plotLegend(legend = c("0.8 - 1.0",
                      "0.6 - 0.8",
                      "0.4 - 0.6",
                      "0.2 - 0.4",
                      "0.0 - 0.2"),
           fill = c("#DD3931", "#EEA741", "#499A53","#98CDED","#262C74"),
           x = 3.6, y = 2.3, width = 0.1, height = 0.4, border = FALSE, 
           fontsize = 6)
plotText(paste0("r2 relative to ", mmp16_sig1_topVar_rsID), x = 4.2, y = 2.2, 
         fontsize = 6, fontfamily = "Helvetica")



# After conditional, signal 2
plotText(label = paste0("FN-f independent eQTL signal 2 conditioned on ", mmp16_sig1_topVar_rsID),
         x = 2.5, y = 3.4, fontfamily = "Helvetica",
         fontsize = 10)


man_condSig2 <- plotManhattan(data = mmp16_fnf_nom_Cond |> filter(signal == 1),
                              params = mmp16_signal_region,
                              fill = colorby("LDgrp",
                                             palette = colorRampPalette(c("#262C74",
                                                                          "#98CDED",
                                                                          "#499A53",
                                                                          "#EEA741",
                                                                          "#DD3931",
                                                                          "grey"))),
                              range = c(0, mmp16_signal_ylim),
                              y = 3.5, height = 1.25,
                              snpHighlights = data.frame(snp = c(mmp16_sig2_topVar_rsID, mmp16_sig1_topVar_rsID),
                                                         pch = c(23, 23), 
                                                         cex = c(0.5, 0.5),
                                                         col = c("#DD3931","#DD3931")))


annoYaxis(plot = man_condSig2, at = seq(0, mmp16_signal_ylim, 2), 
          axisLine = TRUE, fontsize = 8)
plotText(
  label = "-log10(p-value)", x = 0.2, y = 4.15, rot = 90,
  fontsize = 6, just = "center",
  default.units = "inches"
)
plotText(label = mmp16_sig1_topVar_rsID, x = 1.8, y = 4.55,
         fontfamily = "Helvetica", fontsize = 7)
plotText(label = mmp16_sig2_topVar_rsID, x = 3.1, y = 3.95,
         fontfamily = "Helvetica", fontsize = 7)

plotLegend(legend = c("0.8 - 1.0",
                      "0.6 - 0.8",
                      "0.4 - 0.6",
                      "0.2 - 0.4",
                      "0.0 - 0.2"),
           fill = c("#DD3931", "#EEA741", "#499A53","#98CDED","#262C74"),
           x = 3.6, y = 3.7, width = 0.1, height = 0.4, border = FALSE, 
           fontsize = 6)
plotText(paste0("r2 relative to ", mmp16_sig2_topVar_rsID), x = 4.2, y = 3.6, 
         fontsize = 6, fontfamily = "Helvetica")

# Gene track
mmp16_gene <- plotGenes(params = mmp16_signal_region,
                        y = 4.8, height = 0.5, 
                        geneHighlights = data.frame(gene = "MMP16",
                                                    color = "#37a7db"))
annoGenomeLabel(plot = mmp16_gene, y = 5.4, params = mmp16_signal_region,
                fontsize = 8)
dev.off()
