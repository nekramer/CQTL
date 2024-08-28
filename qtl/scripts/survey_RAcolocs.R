library(tidyverse)
library(patchwork)
library(plotgardener)
library(GenomicFeatures)
library(org.Hs.eg.db)
library(grid)
library(RColorBrewer)

yl_gn_bu <- brewer.pal(n = 9, name = "YlGnBu")


## Colocs
ra_colocs <- read_csv("data/RA_PBS_FNF_coloc_results_250kb.csv") |> 
  filter(GWAS_R2 > 0.2 | qtl_R2 > 0.2)

pdf("plots/RA_colocs.pdf", width = 5.5, height = 4)

for (i in 1:nrow(ra_colocs)){
  
  RA_type <- ra_colocs[[i, "RA"]]
  ancestry <- ra_colocs[[i, "ancestry"]]
  if (ancestry == "multi"){
    ancestry <- "multi_ancestry"
  }
  
  GWAS_lead <- ra_colocs[[i, "GWAS_lead"]]
  qtl_variantID <- ra_colocs[[i, "qtl_variantID"]]
  qtl_variant_pos <- as.numeric(unlist(str_split(qtl_variantID, ":"))[2])
  condition <- ra_colocs[[i, "Condition"]]
  if (condition == "PBS"){
    file_cond <- "CTL"
    file_peer <- 20
  } else{
    file_cond <- "FNF"
    file_peer <- 22
  }
  
  # Grab GWAS signal with ld buddies
  gwas_data <- read_csv(paste0("/proj/phanstiel_lab/External/gwas/RA/",
                                     RA_type, "/", ancestry, "/leads/",
                               RA_type, "_", ancestry, "_leads_ld_hg38_final.csv.gz")) |> 
    filter(rsID == GWAS_lead) |> 
    mutate(snp = paste0(gsub("chr", "", chrom), ":", ldbuddy_hg38pos, ":", ldbuddy_NEA, ":", ldbuddy_EA)) |> 
    dplyr::select(variantID_hg38, snp, chrom, ldbuddy_hg38pos, ldbuddy_pval, ldbuddy_R2) |> 
    mutate(ldbuddy_R2 = as.numeric(ldbuddy_R2),
           ldbuddy_hg38pos = as.numeric(ldbuddy_hg38pos)) |> 
    dplyr::rename(pos = ldbuddy_hg38pos,
                  p = ldbuddy_pval)
  
  
  GWAS_variantID <- paste0("chr", unique(gwas_data$variantID_hg38))
  
  # Get summary stats for other variants in region w/o R2 info
  summary_stats <- read_csv(paste0("/proj/phanstiel_lab/External/gwas/RA/",
                                      RA_type, "/", ancestry, "/summary_stats/", 
                                      RA_type, "_", ancestry, "_hg38.csv.gz")) |> 
    filter(chrom == unique(gwas_data$chrom) &
             hg38pos >= qtl_variant_pos - 500000 &
             hg38pos <= qtl_variant_pos + 500000) |>
    filter(!SNP_hg38 %in% gwas_data$snp) |> 
    dplyr::select(SNP_hg38, chrom, hg38pos, Pval) |> 
    mutate(hg38pos = as.numeric(hg38pos)) |> 
    dplyr::rename(snp = SNP_hg38, pos = hg38pos, p = Pval) |> 
    mutate(ldbuddy_R2 = NA)
  
  # Combine for whole signal
  gwas_signal <- gwas_data |> 
    dplyr::select(-variantID_hg38) |> 
    bind_rows(summary_stats) |>  
    filter(!is.na(p)) |> 
    # Create group for LD group
    mutate(LDgrp = factor(cut(ldbuddy_R2, c(0, 0.2, 0.4, 0.6, 0.8, 1)),
                          levels = c("(0,0.2]", "(0.2,0.4]", "(0.4,0.6]",
                                     "(0.6,0.8]", "(0.8,1]", NA)))
  gwas_signal$LDgrp <- addNA(gwas_signal$LDgrp)
  
   
  # eQTL signal
  
  qtl_ld <- read_csv(paste0("/proj/phanstiel_lab/Data/processed/CQTL/eqtl/", 
                            file_cond, "_PEER_k", file_peer, 
                            "_genoPC_cond1Mb_topSignals_rsID_LD.csv")) |>
    filter(gene_id == ra_colocs[[i, "eGene_id"]] & signal == ra_colocs[[i, "signal"]]) |>
    dplyr::select(ld_variantID, R2)
  
  qtl_signal <- 
    read_csv(paste0("/proj/phanstiel_lab/Data/processed/CQTL/eqtl/",
                    file_cond, "_PEER_k", file_peer, "_genoPC_allSignals_nom1Mb_MAFs_",
                    unique(gwas_data$chrom), ".csv")) |> 
    filter(gene_id == ra_colocs[[i, "eGene_id"]] & signal == ra_colocs[[i, "signal"]]) |> 
    mutate(rsID = ifelse(variantID == ra_colocs[[i, "qtl_variantID"]], ra_colocs[[i, "qtl_rsid"]],
                                ifelse(variantID == GWAS_variantID, GWAS_lead, NA))) |> 
    mutate(variant_start = as.numeric(variant_start)) |> 
    dplyr::rename(chrom = variant_chr,
                  pos = variant_start,
                  p = nom_pval,
                  snp = rsID) |> 
    left_join(qtl_ld, by = join_by(variantID == ld_variantID)) |> 
    mutate(LDgrp = cut(R2, c(0, 0.2, 0.4, 0.6, 0.8, 1)))
  qtl_signal$LDgrp <- addNA(qtl_signal$LDgrp)

  
  pageCreate(width = 5.5, height = 4, showGuides = FALSE)
  
  coloc_region <- pgParams(assembly = "hg38",
                           chrom = unique(gwas_data$chrom),
                           chromstart = qtl_variant_pos - 1000000 ,
                           chromend = qtl_variant_pos + 1000000 ,
                           width = 3.9, x = 0.5)
  
  qtl_ylim <- ceiling(max(-1*log10(qtl_signal$p))) + 1.5
  gwas_ylim <- ceiling(max(-1*log10(na.omit(gwas_signal$p)))) + 1.5
  
  # GWAS
  gwas_plot <- plotManhattan(data = gwas_signal |> 
                                     arrange(desc(LDgrp)),
                                   params =  coloc_region,
                                   range = c(0, gwas_ylim),
                                   y = 1.5, height = 1, 
                                   sigline = FALSE, just = c("left", "bottom"),
                                   fill = colorby("LDgrp",
                                                  palette = colorRampPalette(c("#262C74",
                                                                               "#98CDED",
                                                                               "#499A53",
                                                                               "#EEA741",
                                                                               "#DD3931",
                                                                               "grey"))),
                                   snpHighlights = data.frame(snp = c(gsub("chr", "", ra_colocs[[i, "qtl_variantID"]]), 
                                                                      gsub("chr", "", GWAS_variantID)),
                                                              pch = c(24, 23),
                                                              cex = c(0.75, 0.75),
                                                              col = c("black", "black")))
  annoYaxis(plot = gwas_plot, at = seq(0, gwas_ylim, 2), 
            axisLine = TRUE, fontsize = 6)
  
  plotText(
    label = "-log10(p-value)", x = 0.26, y = 1.1, rot = 90,
    fontsize = 6, just = "center",
    default.units = "inches", fontfamily = "Helvetica"
  )
  plotText(label = paste0(ancestry, " ", RA_type, " GWAS"), x = 4.45, y = 0.5, 
           just = c("right", "top"), fontface = "bold",
           fontsize = 8, fontfamily = "Helvetica", lineheight = 0.8)
  
  # eQTL
  qtl_plot <- plotManhattan(data = qtl_signal,
                             params = coloc_region, 
                             range = c(0, qtl_ylim),
                             snpHighlights = data.frame(snp = c(GWAS_lead,
                                                                ra_colocs[[i, "qtl_rsid"]]),
                                                        pch = c(23, 24),
                                                        cex = c(0.75, 0.75),
                                                        col = c("black", "black")),
                             fill = colorby("LDgrp",
                                            palette = colorRampPalette(c("#262C74",
                                                                         "#98CDED",
                                                                         "#499A53",
                                                                         "#EEA741",
                                                                         "#DD3931",
                                                                         "grey"))),
                             y = 1.75, height = 1)
  annoYaxis(plot = qtl_plot, at = seq(0, qtl_ylim, 2), 
            axisLine = TRUE, fontsize = 6)
  plotText(
    label = "-log10(p-value)", x = 0.26, y = 2.45, rot = 90,
    fontsize = 6, just = "center",
    default.units = "inches", fontfamily = "Helvetica"
  )
  plotText(label = paste0(condition, " eQTL"), x = 4.45, y = 1.75, just = "right", fontface = "bold",
           fontsize = 8, fontfamily = "Helvetica")
  plotText(label = paste0("PP4 = ", 
                          signif(ra_colocs[[i, "PP4"]], digits = 3)), 
           x = 4.45, y = 1.9, 
           just = "right", fontsize = 7, fontfamily = "Helvetica")
  
  # Genes
  if (is.na(ra_colocs[[i, "eGene_name"]])){
    ensembl_txdb <- makeTxDbFromEnsembl(organism = "Homo sapiens", release = 111)
    
    gene_plot <- plotGenes(params = coloc_region,
                           chrom = gsub("chr", "", coloc_region$chrom),
                           assembly = assembly(Genome = "ens_hg38",
                                               TxDb = ensembl_txdb,
                                               OrgDb = org.Hs.eg.db, gene.id.column = "ENSEMBL",
                                               display.column = "ENSEMBL"),
                           geneHighlights = data.frame("gene" = ra_colocs[[i, "eGene_id"]],
                                                       "color" = yl_gn_bu[6]),
                           y = 2.9, height = 0.5)
    
  } else {
    gene_plot <- plotGenes(params = coloc_region, 
                           geneHighlights = data.frame("gene" = ra_colocs[[i, "eGene_name"]],
                                                       "color" = yl_gn_bu[6]),
                           y = 2.9, height = 0.5)
  }
  

  annoGenomeLabel(plot = qtl_plot, y = 3.5,
                  x = 0.5, fontsize = 8)
  
  
  
  # LD legend
  grid.points(x = 0.7, y = 3.9, default.units = "native", pch = 24, size = unit(0.45, "char"))
  plotText(label = paste0(ra_colocs[[i, "qtl_rsid"]], " (eQTL index)"),
           fontsize = 7, fontfamily = "Helvetica",
           just = "left", x = 0.8, y = 3.9)
  
  grid.points(x = 2, y = 3.9, default.units = "native", pch = 23, size = unit(0.5, "char"))
  plotText(label = paste0(GWAS_lead," (GWAS index)"),
           fontsize = 7, fontfamily = "Helvetica",
           just = "left", x = 2.1, y = 3.9)
  
  plotLegend(legend = c("0.8 - 1.0",
                        "0.6 - 0.8",
                        "0.4 - 0.6",
                        "0.2 - 0.4",
                        "0.0 - 0.2"),
             fill = c("#DD3931", "#EEA741", "#499A53","#98CDED","#262C74"),
             x = 4.4, y = 0.5, width = 0.1, height = 0.4, border = FALSE, 
             fontsize = 6)
  plotText(label = "r2", x = 5.05, y = 0.4, fontsize = 6)
  
}

dev.off()
