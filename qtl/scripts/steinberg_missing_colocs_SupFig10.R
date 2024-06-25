library(tidyverse)
library(data.table)
library(plotgardener)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(grid)


## NPC1 
# Steinberg: high-grade cartilage eQTL, All OA
# PBS and FN-f eGene
npc1_pbs_egene <- read_csv("data/eqtl/CTL_PEER_k20_genoPC_cond1Mb_topSignals_rsID.csv") |>
  filter(gene_symbol == "NPC1" & signal == 0)
npc1_fnf_egene <- read_csv("data/eqtl/FNF_PEER_k22_genoPC_cond1Mb_topSignals_rsID.csv") |>
  filter(gene_symbol == "NPC1" & signal == 0)
npc1_variant_center <- 0.5*(npc1_pbs_egene[["variant_start"]] + npc1_fnf_egene[["variant_start"]])

npc1_egene_pbs_ld <- fread("data/eqtl/CTL_PEER_k20_genoPC_cond1Mb_topSignals_rsID_LD.csv",
                       data.table = FALSE) |> 
  filter(gene_symbol == "NPC1" & signal == 0) |> 
  dplyr::select(ld_variantID, R2)

npc1_egene_fnf_ld <- fread("data/eqtl/FNF_PEER_k22_genoPC_cond1Mb_topSignals_rsID_LD.csv",
                           data.table = FALSE) |> 
  filter(gene_symbol == "NPC1" & signal == 0) |> 
  dplyr::select(ld_variantID, R2)

# GWAS - npc1 region 
# In Tachmazidou, lead was rs10502437 (chr18:23390742:G:A)

rs10502437_ld <- fread("data/eqtl/qtl_supp/npc1_rs10502437_eur_ld.ld", data.table = FALSE) |> 
  dplyr::select(SNP_B, R2) |> 
  separate_wider_delim(cols = "SNP_B", delim = ":", names = c("chrom", "pos", "A1", "A2"), cols_remove = FALSE) |> 
  mutate(SNP_B_v1 = paste0("chr", SNP_B)) |> 
  # also flip alleles in IDs for attempts at matching with GWAS summary statistics
  mutate(SNP_B_v2 = paste0("chr", chrom, ":", pos, ":", A2, ":", A1))


gwas_npc1 <- fread(paste0("/proj/phanstiel_lab/External/gwas/OA/Boer_2021_hg38_LD/AllOA/summary_stats/AllOA_", 
             npc1_pbs_egene[["gene_chr"]], ".csv"),
      data.table = FALSE) |> 
  filter(hg38pos >= (npc1_variant_center - 500000) & 
           hg38pos <= (npc1_variant_center + 500000)) |> 
  rowwise() |> 
  mutate(snp = paste0("chr", `CHR:hg38POS`, ":", paste(sort(c(EA, NEA)), collapse = ":")),
         chrom = paste0("chr", chrom)) |> 
  dplyr::rename(pos = hg38pos) |> 
  ungroup() |> 
  left_join(rs10502437_ld |> dplyr::select(SNP_B_v1, R2), by = join_by(snp == SNP_B_v1)) |>
  left_join(rs10502437_ld |> dplyr::select(SNP_B_v2, R2), by = join_by(snp == SNP_B_v2)) |> 
  mutate(R2 = case_when(is.na(R2.x) & is.na(R2.y) ~ NA,
                        is.na(R2.x) & !is.na(R2.y) ~ R2.y,
                        !is.na(R2.x) & is.na(R2.y) ~ R2.x)) |> 
  mutate(LDgrp = cut(R2, c(0, 0.2, 0.4, 0.6, 0.8, 1)))

gwas_npc1$LDgrp <- addNA(gwas_npc1$LDgrp)

# eQTL
npc1_signals <- bind_rows(fread(paste0("data/eqtl/qtl_nom/CTL_PEER_k20_genoPC_allSignals_nom1Mb_MAFs_", 
                                       npc1_pbs_egene[["gene_chr"]], ".csv"),
                                data.table = FALSE) |> 
                            filter(gene_id == npc1_pbs_egene[["gene_id"]] & signal ==npc1_pbs_egene[["signal"]]) |> 
                            mutate(Condition = "PBS") |> 
                            left_join(npc1_egene_pbs_ld, by = join_by(variantID == ld_variantID)),
                          fread(paste0("data/eqtl/qtl_nom/FNF_PEER_k22_genoPC_allSignals_nom1Mb_MAFs_", 
                                       npc1_pbs_egene[["gene_chr"]], ".csv"),
                                data.table = FALSE) |> 
                            filter(gene_id == npc1_pbs_egene[["gene_id"]] & signal == npc1_pbs_egene[["signal"]]) |> 
                            mutate(Condition = "FN-f") |> 
                            left_join(npc1_egene_fnf_ld, by = join_by(variantID == ld_variantID))) |>
  dplyr::rename(chrom = variant_chr,
                pos = variant_start,
                p = nom_pval,
                snp = variantID) |> 
  mutate(LDgrp = cut(R2, c(0, 0.2, 0.4, 0.6, 0.8, 1)))
npc1_signals$LDgrp <- addNA(npc1_signals$LDgrp)


## FAM53A (ENSG00000174137)
# Steinberg: low-grade cartilage eQTL, KneeHip OA
# not an eGene in either PBS or FN-f
fam53a_pbs_variant <- read_csv("data/eqtl/qtl_supp/CTL_PEER_k20_genoPC_perm1Mb_FDR.csv") |> 
  filter(gene_symbol == "FAM53A")
fam53a_fnf_variant <- read_csv("data/eqtl/qtl_supp/FNF_PEER_k22_genoPC_perm1Mb_FDR.csv") |> 
  filter(gene_symbol == "FAM53A")


# GWAS - fam53a region, rsID rs1530586
#In Tachmazidou, lead was rs11732213 (chr4:1702517:T:C)
# Boer reported signal to be associated with TACC3 in TJR?
fam53a_variant_center <- 1702517
rs11732213_ld <- fread("data/eqtl/qtl_supp/fam53a_rs11732213_eur_ld.ld", data.table = FALSE) |> 
  dplyr::select(SNP_B, R2) |> 
  separate_wider_delim(cols = "SNP_B", delim = ":", names = c("chrom", "pos", "A1", "A2"), cols_remove = FALSE) |> 
  mutate(SNP_B_v1 = paste0("chr", SNP_B)) |> 
  # also flip alleles in IDs for attempts at matching with GWAS summary statistics
  mutate(SNP_B_v2 = paste0("chr", chrom, ":", pos, ":", A2, ":", A1))

gwas_fam53a <- fread("/proj/phanstiel_lab/External/gwas/OA/Boer_2021_hg38_LD/KneeHipOA/summary_stats/KneeHipOA_chr4.csv",
                   data.table = FALSE) |> 
  filter(hg38pos >= fam53a_variant_center - 500000 &
           hg38pos <= fam53a_variant_center + 500000) |> 
  rowwise() |> 
  mutate(snp = paste0("chr", `CHR:hg38POS`, ":", paste(sort(c(EA, NEA)), collapse = ":")),
         chrom = paste0("chr", chrom)) |> 
  dplyr::rename(pos = hg38pos) |> 
  ungroup() |> 
  left_join(rs11732213_ld |> dplyr::select(SNP_B_v1, R2), by = join_by(snp == SNP_B_v1)) |>
  left_join(rs11732213_ld |> dplyr::select(SNP_B_v2, R2), by = join_by(snp == SNP_B_v2)) |> 
  mutate(R2 = case_when(is.na(R2.x) & is.na(R2.y) ~ NA,
                        is.na(R2.x) & !is.na(R2.y) ~ R2.y,
                        !is.na(R2.x) & is.na(R2.y) ~ R2.x)) |> 
  mutate(LDgrp = cut(R2, c(0, 0.2, 0.4, 0.6, 0.8, 1)))

gwas_fam53a$LDgrp <- addNA(gwas_fam53a$LDgrp)

# eQTL 
fam53a_signals <- bind_rows(fread("data/eqtl/qtl_nom/CTL_PEER_k20_genoPC_allSignals_nom1Mb_MAFs_chr4.csv",
                                data.table = FALSE) |> 
                              filter(gene_id == "ENSG00000174137" & signal == 0) |> 
                              mutate(Condition = "PBS") |> 
                              left_join(rs11732213_ld |> dplyr::select(SNP_B_v1, R2), by = join_by(variantID == SNP_B_v1)) |> 
                              left_join(rs11732213_ld |> dplyr::select(SNP_B_v2, R2), by = join_by(variantID == SNP_B_v2)) |> 
                              mutate(R2 = case_when(is.na(R2.x) & is.na(R2.y) ~ NA,
                                                    is.na(R2.x) & !is.na(R2.y) ~ R2.y,
                                                    !is.na(R2.x) & is.na(R2.y) ~ R2.x)),
                          fread("data/eqtl/qtl_nom/FNF_PEER_k22_genoPC_allSignals_nom1Mb_MAFs_chr4.csv",
                                data.table = FALSE) |> 
                            filter(gene_id == "ENSG00000174137" & signal == 0) |> 
                            mutate(Condition = "FN-f") |> 
                            left_join(rs11732213_ld |> dplyr::select(SNP_B_v1, R2), by = join_by(variantID == SNP_B_v1)) |> 
                            left_join(rs11732213_ld |> dplyr::select(SNP_B_v2, R2), by = join_by(variantID == SNP_B_v2)) |> 
                            mutate(R2 = case_when(is.na(R2.x) & is.na(R2.y) ~ NA,
                                                  is.na(R2.x) & !is.na(R2.y) ~ R2.y,
                                                  !is.na(R2.x) & is.na(R2.y) ~ R2.x))) |>
  dplyr::rename(chrom = variant_chr,
                pos = variant_start,
                p = nom_pval,
                snp = variantID) |> 
  mutate(LDgrp = cut(R2, c(0, 0.2, 0.4, 0.6, 0.8, 1)))

fam53a_signals$LDgrp <- addNA(fam53a_signals$LDgrp)


# Plot 
#pdf(file = "plots/colocs_Fig5_supp/SupFig10.pdf", width = 8.6, height = 4.75)
pageCreate(width = 8.6, height = 4.75, showGuides = TRUE)

## A - NPC1 (TMEM241)

plotText("A", x = 0, y = 0.1641, just = "left", fontfamily = "Helvetica",
         fontsize = 9, fontface = "bold")

npc1_region <- pgParams(assembly = "hg38",
                        chrom = npc1_pbs_egene[["gene_chr"]],
                        chromstart = as.integer(npc1_variant_center - 500000),
                        chromend = as.integer(npc1_variant_center + 500000),
                        x = 0.6, width = 3.25)


qtl_npc1_ylim <- ceiling(max(-1*log10(npc1_signals$p))) + 1.5

gwas_npc1_plot <- plotManhattan(data = gwas_npc1 |> arrange(desc(LDgrp)),
                           params = npc1_region,
                           range = c(0, qtl_npc1_ylim),
                           y = 0.15, height = 1.15, 
                           sigLine = TRUE, lty = 2,
                           fill = colorby("LDgrp",
                                          palette = colorRampPalette(c("#262C74",
                                                                       "#98CDED",
                                                                       "#499A53",
                                                                       "#EEA741",
                                                                       "#DD3931",
                                                                       "grey"))),
                           snpHighlights = data.frame(snp = c(npc1_pbs_egene[["variantID"]], 
                                                              npc1_fnf_egene[["variantID"]], 
                                                              "chr18:23390742:A:G"),
                                                      pch = c(24, 25 ,23),
                                                      cex = c(0.75, 0.75, 0.75),
                                                      col = c("white", "white", "black")))

annoYaxis(plot = gwas_npc1_plot, at = seq(0, qtl_npc1_ylim, 2), 
          axisLine = TRUE, fontsize = 7)
plotText(
  label = "-log10(p-value)", x = 0.2525, y = 0.75, rot = 90,
  fontsize = 7, just = "center",
  default.units = "inches", fontfamily = "Helvetica"
)

plotText(label = "All OA GWAS", x = 3.85, y = 0.2, just = "right", 
         fontsize = 10, fontfamily = "Helvetica")
pbs_npc1 <- plotManhattan(data = npc1_signals |> filter(Condition == "PBS"),
                          params = npc1_region, 
                          range = c(0, qtl_npc1_ylim),
                          snpHighlights = data.frame(snp = c(npc1_pbs_egene[["variantID"]], npc1_fnf_egene[["variantID"]],
                                                             "chr18:23390742:G:A"),
                                                     pch = c(24, 25, 23),
                                                     cex = c(0.75, 0.75, 0.75),
                                                     col = c("black", "black", "white")),
                          fill = colorby("LDgrp",
                                         palette = colorRampPalette(c("#262C74",
                                                                      "#98CDED",
                                                                      "#499A53",
                                                                      "#EEA741",
                                                                      "#DD3931",
                                                                      "grey"))),
                          y = 1.3, height = 1.15)

annoYaxis(plot = pbs_npc1, at = seq(0, qtl_npc1_ylim, 2), 
          axisLine = TRUE, fontsize = 7)
plotText(
  label = "-log10(p-value)", x = 0.2525, y = 2, rot = 90,
  fontsize = 7, just = "center",
  default.units = "inches", fontfamily = "Helvetica"
)
plotText(label = "PBS eQTL", x = 3.85, y = 1.45, just = "right", fontface = "bold",
         fontsize = 10, fontfamily = "Helvetica", fontcolor = "#6B98C9")

fnf_npc1 <- plotManhattan(data = npc1_signals |> filter(Condition == "FN-f"),
                          params = npc1_region, 
                          range = c(0, qtl_npc1_ylim),
                          snpHighlights = data.frame(snp = c(npc1_pbs_egene[["variantID"]], npc1_fnf_egene[["variantID"]],
                                                             "chr18:23390742:G:A"),
                                                     pch = c(24, 25, 23),
                                                     cex = c(0.75, 0.75, 0.75),
                                                     col = c("black", "black", "white")),
                          fill = colorby("LDgrp",
                                         palette = colorRampPalette(c("#262C74",
                                                                      "#98CDED",
                                                                      "#499A53",
                                                                      "#EEA741",
                                                                      "#DD3931",
                                                                      "grey"))),
                          y = 2.45, height = 1.15)

annoYaxis(plot = fnf_npc1, at = seq(0, qtl_npc1_ylim, 2), 
          axisLine = TRUE, fontsize = 7)
plotText(
  label = "-log10(p-value)", x = 0.2525, y = 3.1, rot = 90,
  fontsize = 7, just = "center",
  default.units = "inches", fontfamily = "Helvetica"
)
plotText(label = "FN-f eQTL", x = 3.85, y = 2.6, just = "right", fontface = "bold",
         fontsize = 10, fontfamily = "Helvetica", fontcolor = "#D6912B")

npc1_genes <- plotGenes(params = npc1_region, y = 3.65,
                        height = 0.5, geneHighlights = data.frame("gene" = c("NPC1", "TMEM241"),
                                                                  "color" = c("#37a7db", "lightsteelblue2")))
annoGenomeLabel(plot = npc1_genes, params = npc1_region, y = 4.15, fontsize = 8,
                commas = TRUE)
grid.points(x = 0.6, y = 4.4, default.units = "native", pch = 23, size = unit(0.75, "char"))
plotText(label = "rs10502437 (Tachmazidou et al. GWAS index)",
         x = 0.7, y = 4.4, fontfamily = "Helvetica", fontsize = 6, just = "left")

grid.points(x = 2.6, y = 4.4, default.units = "native", pch = 24, size = unit(0.65, "char"))
plotText(label = "rs8083301 (PBS eQTL index)",
         x = 2.7, y = 4.4, fontfamily = "Helvetica", fontsize = 6, just = "left")

grid.points(x = 1.6, y = 4.55, default.units = "native", pch = 25, size = unit(0.65, "char"))
plotText(label = "rs6507716 (FN-f eQTL index)",
         x = 1.75, y = 4.55, fontfamily = "Helvetica", fontsize = 6, just = "left")

## B - FAM53A (SLBP)

plotText("B", x = 4, y = 0.1, just = c("left", "top"), fontfamily = "Helvetica",
         fontsize = 14, fontface = "bold")

fam53a_plot_region <- pgParams(assembly = "hg38",
                               chrom = fam53a_pbs_variant[["gene_chr"]],
                               chromstart = as.integer(fam53a_variant_center - 500000),
                               chromend = as.integer(fam53a_variant_center + 500000),
                        x = 4.5, width = 3.25)


gwas_fam53a_ylim <- ceiling(max(-1*log10(gwas_fam53a$p))) + 1.5

gwas_fam53a_plot <- plotManhattan(data = gwas_fam53a |> arrange(desc(LDgrp)),
                                params = fam53a_plot_region,
                                range = c(0, gwas_fam53a_ylim),
                                y = 0.15, height = 1.15, 
                                sigLine = TRUE, lty = 2,
                                fill = colorby("LDgrp",
                                               palette = colorRampPalette(c("#262C74",
                                                                            "#98CDED",
                                                                            "#499A53",
                                                                            "#EEA741",
                                                                            "#DD3931",
                                                                            "grey"))),
                                snpHighlights = data.frame(snp = c("chr4:1261672:A:G", fam53a_fnf_variant[["variantID"]],
                                                                   "chr4:1702517:C:T", "chr4:1759200:C:T"),
                                                           pch = c(24, 25, 23, 22),
                                                           cex = c(0.75, 0.75, 0.75, 0.75),
                                                           col = c("black", "white", "black", "black")))

annoYaxis(plot = gwas_fam53a_plot, at = seq(0, gwas_fam53a_ylim, 2), 
          axisLine = TRUE, fontsize = 7)
plotText(
  label = "-log10(p-value)", x = 4.125, y = 0.75, rot = 90,
  fontsize = 7, just = "center",
  default.units = "inches", fontfamily = "Helvetica"
)
plotText(label = "KneeHip OA GWAS", x = 7.75, y = 0.2, just = "right", fontface = "bold",
         fontsize = 10, fontfamily = "Helvetica")
pbs_fam53a <- plotManhattan(data = fam53a_signals |> filter(Condition == "PBS"),
                          params = fam53a_plot_region, 
                          range = c(0, gwas_fam53a_ylim),
                          y = 1.3, height = 1.15, 
                          fill = colorby("LDgrp",
                                                palette = colorRampPalette(c("#262C74",
                                                                             "#98CDED",
                                                                             "#499A53",
                                                                             "#EEA741",
                                                                             "#DD3931",
                                                                             "grey"))),
                          snpHighlights = data.frame(snp = c("chr4:1702517:T:C", "chr4:1759200:C:T"),
                                                     pch = c(23, 22),
                                                     cex = c(0.75, 0.75), 
                                                     col = c("black", "black")))

annoYaxis(plot = pbs_fam53a, at = seq(0, gwas_fam53a_ylim, 2), 
          axisLine = TRUE, fontsize = 7)
plotText(
  label = "-log10(p-value)", x = 4.125, y = 2, rot = 90,
  fontsize = 7, just = "center",
  default.units = "inches", fontfamily = "Helvetica"
)
plotText(label = "PBS eQTL", x = 7.75, y = 1.45, just = "right", fontface = "bold",
         fontsize = 10, fontfamily = "Helvetica", fontcolor = "#6B98C9")

fnf_fam53a <- plotManhattan(data = fam53a_signals |> filter(Condition == "FN-f"),
                          params = fam53a_plot_region, 
                          range = c(0, gwas_fam53a_ylim),
                          y = 2.45, height = 1.15,
                          fill = colorby("LDgrp",
                                         palette = colorRampPalette(c("#262C74",
                                                                      "#98CDED",
                                                                      "#499A53",
                                                                      "#EEA741",
                                                                      "#DD3931",
                                                                      "grey"))),
                          snpHighlights = data.frame(snp = c("chr4:1702517:T:C", "chr4:1759200:C:T"),
                                                                      pch = c(23, 22),
                                                                      cex = c(0.75, 0.75),
                                                                      col = c("black", "black")))

annoYaxis(plot = fnf_fam53a, at = seq(0, gwas_fam53a_ylim, 2), 
          axisLine = TRUE, fontsize = 7)
plotText(
  label = "-log10(p-value)", x = 4.125, y = 3.1, rot = 90,
  fontsize = 7, just = "center",
  default.units = "inches", fontfamily = "Helvetica"
)
plotText(label = "FN-f eQTL", x = 7.75, y = 2.6, just = "right", fontface = "bold",
         fontsize = 10, fontfamily = "Helvetica", fontcolor = "#D6912B")

fam53a_genes <- plotGenes(params = fam53a_plot_region, y = 3.65,
                        height = 0.5, geneHighlights = data.frame("gene" = c("FAM53A", "SLBP", "TACC3"),
                                                                  "color" = c("#37a7db","lightsteelblue2", "grey")),
                        fontsize = 4)
annoGenomeLabel(plot = fam53a_genes, params = fam53a_plot_region, y = 4.15, fontsize = 8,
                commas = TRUE)

grid.points(x = 4.5, y = 4.4, default.units = "native", pch = 23, size = unit(0.75, "char"))
plotText(label = "rs11732213 (Tachmazidou et al. GWAS index)",
         x = 4.6, y = 4.4, fontfamily = "Helvetica", fontsize = 6, just = "left")

grid.points(x = 6.375, y = 4.4, default.units = "native", pch = 22, size = unit(0.65, "char"))
plotText(label = "rs1530586 (Boer et al. GWAS index)",
         x = 6.45, y = 4.4, fontfamily = "Helvetica", fontsize = 6, just = "left")

legend <- plotLegend(legend = c("0.8 - 1.0",
                      "0.6 - 0.8",
                      "0.4 - 0.6",
                      "0.2 - 0.4",
                      "0.0 - 0.2"),
           fill = c("#DD3931", "#EEA741", "#499A53","#98CDED","#262C74"),
           x = 7.7, y = 0.25, width = 0.1, height = 0.4, border = FALSE, 
           fontsize = 6)
grid.draw(legend$grobs)
plotText(label = "r2", x = 8.35, y = 0.15, fontsize = 6, fontfamily = "Helvetica")

dev.off()

