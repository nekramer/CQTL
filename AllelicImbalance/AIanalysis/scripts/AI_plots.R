library(VennDiagram)
library(DESeq2)
library(RColorBrewer)
library(data.table)
library(tidyr)
library(ggplot2)
source("scripts/utils.R")
source("scripts/plotAIdata.R")
# Read in and format CTL/FNF AI rsids/genes ------------------------------------

CTL <- reformatRSIDdata("data/2022-08-24_AIsigCTL_rsids.csv")
FNF <- reformatRSIDdata("data/2022-08-24_AIsigFNF_rsids.csv")

# DESeq results --------------------------------------------------------
load("data/2022-08-23_AIsigCTL.rda")
load("data/2022-08-23_AIsigFNF.rda")
load("data/dds.rda")

# Reformat and join with rsid/gene info
sigCTL_df <- reformatDESeqdata(sigCTL, CTL)
sigFNF_df <- reformatDESeqdata(sigFNF, FNF)


# CTL/FNF/both comparisons -----------------------------------------------------

CTL_snps_only <- sigCTL_df[which(!sigCTL_df$rsid %in% sigFNF_df$rsid),]
FNF_snps_only <- sigFNF_df[which(!sigFNF_df$rsid %in% sigCTL_df$rsid),]
both_snps_only <- sigCTL_df[which(sigCTL_df$rsid %in% sigFNF_df$rsid),]


CTL_genes_only <- sigCTL_df[which(!sigCTL_df$gene_symbol %in% sigFNF_df$gene_symbol),]
FNF_genes_only <- sigFNF_df[which(!sigFNF_df$gene_symbol %in% sigFNF_df$gene_symbol),]
both_genes_only <- sigCTL_df[which(sigCTL_df$gene_symbol %in% sigFNF_df$gene_symbol),]

# Venn diagram of # SNPs --------------------------------------------------
venn.diagram(list(sigCTL_df$rsid, sigFNF_df$rsid), filename = "plots/snpVenn.tiff",
             category.names = c("Control", "FN-f"),
             output = TRUE,
             lwd = 0, 
             fill = c("#7bc5ee", "#fc7971"),
             label.col = "grey35",
             cat.col = c("#7bc5ee", "#fc7971"),
             cat.pos = c(360, 360), 
             cex = 2,
             cat.cex = 2.1,
             rotation.degree = 180)

venn.diagram(list(sigCTL_df$gene_symbol, sigFNF_df$gene_symbol), filename = "plots/geneVenn.tiff",
             category.names = c("Control", "FN-f"),
             output = TRUE,
             lwd = 0, 
             fill = c("#7bc5ee", "#fc7971"),
             label.col = "grey35",
             cat.col = c("#7bc5ee", "#fc7971"),
             cat.pos = c(360, 360), 
             cex = 2,
             cat.cex = 2.1,
             rotation.degree = 180)

# AI slope plots -------------------------------------------------------------
plotSNPs <- function(data, dds, save_category){
  snp_name <- paste0(data["chr"], ":", data["pos"], ":",
                     data["ref"], ":", data["alt"])
  snp_name <- gsub(" ", "", snp_name)
  rsid <- data["rsid"]
  gene <- data["gene_symbol"]
  AIplot <- plotAIdata(snp = snp_name, gene = gene, dds = dds, rsid = rsid)
  ggsave(AIplot, filename = paste0("plots/", save_category, "/", gene, "_AI.pdf"))
}

apply(both_snps_only, 1, plotSNPs, dds = dds, save_category= "both")
apply(CTL_snps_only, 1, plotSNPs, dds = dds, save_category = "CTL")
apply(FNF_snps_only, 1, plotSNPs, dds = dds, save_category = "FNF")

# allele count fractions --------------------------------------------------

plotAltFractions <- function(data, dds, save_category){
  snp_name <- paste0(data["chr"], ":", data["pos"], ":",
                     data["ref"], ":", data["alt"])
  snp_name <- gsub(" ", "", snp_name)
  rsid <- data["rsid"]
  gene <- data["gene_symbol"]
  altplot <- plotDonorAltFraction(snp = snp_name, gene = gene, dds = dds, rsid = rsid)
  if (!is.null(altplot)){
    ggsave(altplot, filename = paste0("plots/", save_category, "/", gene, "_altCounts.pdf"))
  }
  
}

apply(both_snps_only, 1, plotAltFractions, dds = dds, save_category = "both")
apply(CTL_snps_only, 1, plotAltFractions, dds = dds, save_category = "CTL")
apply(FNF_snps_only, 1, plotAltFractions, dds = dds, save_category = "FNF")
  

colData <- fread("../processing/output/AI/colData.csv", data.table = FALSE)
colnames <- paste0(colData$donor, "_", colData$condition, "_", colData$allele)
alleleCountsMatrix <- fread("../processing/output/AI/alleleCountsMatrix.csv", data.table = FALSE,
                            colClasses = c("character", rep("numeric", nrow(colData))),
                            skip = 4) %>% separate(V1, c("chr", "pos", "ref", "alt"))
colnames(alleleCountsMatrix)[5:ncol(alleleCountsMatrix)] <- colnames











