library(data.table)
library(VennDiagram)
library(tidyverse)
library(ggrepel)

# Join with LD buddies for each SNP ------------------------------------------

# CTL_LDbuddies <- read_csv(file = "data/CTL_LDbuddies.csv")
# FNF_LDbuddies <- read_csv(file = "data/FNF_LDbuddies.csv")
# 
# sigCTL_df_alleleCounts_ldbuddies <- left_join(sigCTL_df_alleleCounts, CTL_LDbuddies, by = "variantID")
# sigFNF_df_alleleCounts_ldbuddies <- left_join(sigFNF_df_alleleCounts, FNF_LDbuddies, by = "variantID")


AI_results <- fread("data/AIresults_final.csv", data.table = FALSE)

# Functions ---------------------------------------------------------------

# Function to plot a SNP's AI with a slope plot between alleles
makeAISlopePlot <- function(snp, data, sigGroup){
  
  # Subset data for SNP and filter for heterozygous donors
  data_subset <- data %>% filter(rsid == snp) %>%
    group_by(donor, condition) %>% 
    filter(sum(count) >= 10) %>% 
    filter(all(count >= 2)) %>%
    ungroup()
  
  geneName <- data_subset %>% pull(gene_symbol) %>% unique()
  geneName <- str_c(geneName, collapse = "_")
  # Set allele factors
  data_subset$allele <- factor(data_subset$allele, levels = c("ref", "alt"))
  refAllele <- data_subset %>% filter(allele == "ref") %>% pull(ref) %>% unique()
  altAllele <- data_subset %>% filter(allele == "alt") %>% pull(alt) %>% unique()
  
  plot <- ggplot(data = data_subset, aes(allele, count, color = condition)) + 
    geom_text_repel(data = unique(subset(data_subset, allele == "alt")), 
                    aes(label = donor), x ="alt", 
                    direction = "y", hjust = -0.35,
                    min.segment.length = 0.1) +
    geom_text_repel(data = unique(subset(data_subset, allele == "ref")), 
                    aes(label = donor), x ="ref", 
                    direction = "y", hjust = 1.5,
                    min.segment.length = 0.1) +
    scale_color_manual(values = c("#48A9DF", "#EF734A")) +
    scale_x_discrete(labels = c(refAllele, altAllele)) +
    theme_minimal() +
    geom_hline(aes(yintercept = 0), color = "grey") +
    geom_line(aes(group = donor), color = "grey25", lwd = 0.2) +
    geom_point(size = 2) + 
    facet_wrap(~condition, strip.position = "bottom", 
               labeller = as_labeller(c("CTL" = "Control",
                                        "FNF" = "FN-f"))) +
    ylab("normalized RNA-seq counts") +
    theme(panel.grid = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 14),
          axis.text.x.bottom = element_text(vjust = 15, size = 14),
          strip.text.x = element_text(vjust = 2, size = 16),
          legend.position = "none",
          axis.text.y = element_text(size = 12),
          axis.ticks.y = element_line(color = "grey"),
          plot.title = element_text(face = "bold", size = 16),
          plot.subtitle = element_text(size = 14)) +
    labs(title = snp,
         subtitle = geneName)
  
  ggsave(plot = plot, filename = paste0("plots/slopePlots/", sigGroup, "/", snp, "_" , geneName, "_AIslopeplot.pdf"),
         units = "in",
         width = 8,
         height = 8)
  
}

# Function to plot a SNP's AI with alt allele fraction per donor
makeAltFractionPlot <- function(snp, data, sigGroup){
  
  # Subset data for SNP, calculate alt fractions, and subset for heterozygous donors
  data_subset <- data %>% filter(rsid == snp) %>%
    group_by(donor, condition) %>% 
    mutate(total = sum(count)) %>%
    filter(total >= 10) %>% 
    filter(all(count >= 2)) %>%
    ungroup() %>%
    filter(allele == "alt") %>%
    mutate(alt_fraction = count/total)
  
  geneName <- data_subset %>% pull(gene_symbol) %>% unique()
  geneName <- str_c(geneName, collapse = "_")
  
  
  if (sigGroup == "CTL"){
    shapes <- c(17, 16)
  } else if (sigGroup == "FNF"){
    shapes <- c(16, 17)
  } else if (sigGroup == "both"){
    shapes <- c(17, 17)
  }
  
 
  plot <- ggplot(data = data_subset, aes(x = donor, y = alt_fraction, color = condition)) + 
    scale_color_manual(values = c("#48A9DF", "#EF734A"), labels = c("Control", "FN-f"), guide = "none") +
    geom_hline(aes(yintercept = 0.5), color = "grey", linetype = 2) +
    geom_point(size = 3, aes(shape = condition)) + 
    scale_shape_manual(values = shapes, guide = "none") +
    theme_minimal() +
    scale_y_continuous(breaks = c(0, 0.25, 0.50, 0.75, 1), 
                       limits = c(0,1),
                       labels = c(0, 0.25, 0.50, 0.75, 1))+
    ylab("alt allele fraction") +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 45),
          axis.title.y = element_text(size = 14),
          axis.text.y = element_text(size = 12),
          legend.title = element_blank(), 
          plot.title = element_text(face = "bold", size = 16),
          plot.subtitle = element_text(size = 14)) +
    labs(title = snp, subtitle = geneName) 

  
  ggsave(plot = plot, filename = paste0("plots/altFractionPlots/", sigGroup, "/", snp, "_" , geneName, "_AIaltFraction.pdf"),
         units = "in",
         width = 8,
         height = 8)
  
}

# Significant AI ----------------------------------------------------------

# Significant if padj < 0.05 and abs(log2FoldChange) > log2(1.1)
# SNPs
sig_CTL_snps <- AI_results %>%
  filter(CTL_padj < 0.05 & abs(CTL_log2FoldChange) > log2(1.1)) %>%
  distinct(rsid, .keep_all = TRUE)

sig_FNF_snps <- AI_results %>%
  filter(FNF_padj < 0.05 & abs(FNF_log2FoldChange) > log2(1.1)) %>%
  distinct(rsid, .keep_all = TRUE)

# sig_CTL_snps <- AI_results %>%
#   filter(CTL_padj < 0.05) %>%
#   distinct(rsid, .keep_all = TRUE)
# 
# sig_FNF_snps <- AI_results %>%
#   filter(FNF_padj < 0.05) %>%
#   distinct(rsid, .keep_all = TRUE)


CTL_snps_only <- sig_CTL_snps[which(!sig_CTL_snps$rsid %in% sig_FNF_snps$rsid),]
FNF_snps_only <- sig_FNF_snps[which(!sig_FNF_snps$rsid %in% sig_CTL_snps$rsid),]
both_snps_only <- sig_CTL_snps[which(sig_CTL_snps$rsid %in% sig_FNF_snps$rsid),]

# Genes
sig_CTL_genes <- AI_results %>%
  filter(CTL_padj < 0.05 & abs(CTL_log2FoldChange) > log2(1.1)) %>%
  distinct(gene_symbol, .keep_all = TRUE)

sig_FNF_genes <- AI_results %>%
  filter(FNF_padj < 0.05 & abs(FNF_log2FoldChange) > log2(1.1)) %>%
  distinct(gene_symbol, .keep_all = TRUE)


CTL_genes_only <- sig_CTL_genes[which(!sig_CTL_genes$gene_symbol %in% sig_FNF_genes$gene_symbol),]
FNF_genes_only <- sig_FNF_genes[which(!sig_FNF_genes$gene_symbol %in% sig_CTL_genes$gene_symbol),]
both_genes_only <- sig_CTL_genes[which(sig_CTL_genes$gene_symbol %in% sig_FNF_genes$gene_symbol),]

# Venn Diagrams ------------------------------------

# SNPs
venn.diagram(list(sig_CTL_snps$rsid, sig_FNF_snps$rsid), filename = "plots/sigAI_snps_VennDiagram.tiff",
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

# Genes
venn.diagram(list(sig_CTL_genes$gene_symbol, sig_FNF_genes$gene_symbol), filename = "plots/sigAI_genes_VennDiagram.tiff",
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


# Alt allele fraction distributions ---------------------------------------


# Subset data for significance categories ---------------------------------

# sig SNPs only found in CTL
CTL_sigonly_data <- AI_results %>%
  filter(rsid %in% CTL_snps_only$rsid)


# sig SNPs only found in FNF
FNF_sigonly_data <- AI_results %>%
  filter(rsid %in% FNF_snps_only$rsid)


# sig SNPs found in CTL and FNF
both_sig_data <- AI_results %>%
  filter(rsid %in% both_snps_only$rsid)

# Slope plots ------------------------------------------------------------

lapply(CTL_snps_only$rsid, makeAISlopePlot, data = CTL_sigonly_data, sigGroup = "CTL")
lapply(FNF_snps_only$rsid, makeAISlopePlot, data = FNF_sigonly_data, sigGroup = "FNF")
lapply(both_snps_only$rsid, makeAISlopePlot, data = both_sig_data, sigGroup = "both")


# Alt allele fraction plots -----------------------------------------------

lapply(CTL_snps_only$rsid, makeAltFractionPlot, data = CTL_sigonly_data, sigGroup = "CTL")
lapply(FNF_snps_only$rsid, makeAltFractionPlot, data = FNF_sigonly_data, sigGroup = "FNF")
lapply(both_snps_only$rsid, makeAltFractionPlot, data = both_sig_data, sigGroup = "both")

