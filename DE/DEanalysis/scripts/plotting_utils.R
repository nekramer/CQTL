library(ggplot2)
library(DESeq2)
library(tidyverse)
library(grid)
library(colorspace)



conditionColors <- c("CTL" = "#B8B8B8", "FNF" = "#4A4A4A")
ageColors <- sequential_hcl(n = 6, palette = "Mint")
raceColors <- c("#C74A53", "#EFB06E", "#5BAD58", "#177F97", "#704D9E")
sexColors <- c("F" = "#DD8492", "M" = "#4788BA")

heatmapColors <- colorRampPalette(c("#4FAFE1", "black", "#F5D348"))(7)
log2fcColors <- c("+" = "#FBBE67", "-" = "#78A1Cd")


ageClusterColors <- c("1" = "#009E8E", "2" = "#AC7F21")


plotConditionExpression <- function(gene, dds){
  
  geneID <- gene[["gene_id"]]
  geneSymbol <- gene[["symbol"]]
  log2FC <- gene[["log2FoldChange"]]
  
  gene_counts <- plotCounts(dds, gene = geneID, 
                            intgroup = "Condition", 
                            returnData = TRUE) %>%
    rownames_to_column(var = "Sample") %>%
    separate("Sample", 
             into = c(NA, "donor", NA, NA, NA, NA),
             sep = "_", 
             remove = FALSE)
  
  ggplot(gene_counts, aes(Condition, count)) +
    geom_jitter(position = position_jitter(width = 0.1, height = 0.1)) +
    scale_y_continuous(trans = "log") +
    scale_x_discrete(labels = c("Control", "FN-f")) +
    ylab("normalized RNA-seq counts") +
    theme_light() +
    labs(title = geneSymbol, 
         subtitle = paste0("log2FC = ", 
                           round(log2FC, 
                                 digits = 2))) +
    theme(axis.title.x = element_blank(), 
          plot.title = element_text(hjust = 0.5))
  ggsave(filename = paste0("plots/", geneSymbol, "_Conditioncounts.pdf"), 
         units = "in",
         width = 5, height = 6)
  
}

theme_custom_general <- function(){
  
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(), 
          axis.ticks.x = element_blank(),
          panel.background = element_blank(), 
          axis.line.y = element_line(color = "grey25", 
                                     linewidth = 0.25),
          axis.ticks.length.y = unit(-0.1, "cm"),
          axis.ticks.y = element_line(color = "grey25",
                                      linewidth = 0.25))
  
}


theme_custom_scatterplot <- function(){
  
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(color = "grey25", 
                                   linewidth = 0.25),
        axis.ticks.length = unit(-0.1, "cm"),
        axis.ticks = element_line(color = "grey25",
                                    linewidth = 0.25),
        text = element_text(family = "Helvetica"),
        legend.key = element_blank(),
        legend.background = element_blank())
  
}

reduceGO <- function(goterms, category, ont = "BP", threshold = 0.7){
  
  # Calculate similarity matrix for GO terms based on ontology
  simMatrix <- calculateSimMatrix(goterms$TermID,
                                  orgdb = "org.Hs.eg.db",
                                  ont = ont,
                                  method = "Rel")
  
  # Create named vector of scores
  # Higher is better -> -log10 transforming p-values
  scores <- setNames(-log10(goterms$pval), goterms$TermID)
  
  # Group GO terms based on similarity threshold
  reducedTerms <- reduceSimMatrix(simMatrix,
                                  scores,
                                  threshold = threshold,
                                  orgdb = "org.Hs.eg.db") %>%
    distinct(parentTerm, .keep_all = TRUE) %>%
    dplyr::rename(`-log10pval` = score) %>%
    mutate("category" = category)
  
  return(reducedTerms)
}



venn_font <- function(p, font){
  
  grep_grob <- function(gt, lab){
    which(sapply(gt, function(x) grepl(lab, x$name)))
  }
  
  p2 <- ggplot_gtable(ggplot_build(p))
  mygrobs <- p2$grobs
  # Break down grobs into panel, venn object, and text pieces
  panel_grob <- mygrobs[[grep_grob(mygrobs, "panel")]]
  venn_grob <- panel_grob$children[[grep_grob(panel_grob$children, "venn")]]
  text_grobs <- venn_grob$children[grep_grob(venn_grob$children, "text")]
  # Make both new font family
  text_grobs <- do.call(grid::gList, 
                        lapply(text_grobs, 
                               function(x) {x$gp$fontfamily <- font; 
                               x}))
  # Make titles bold
  text_grobs[[1]]$gp$fontface <- "bold"
  
  # Add grobs back
  venn_grob$children[grep_grob(venn_grob$children, "text")] <- text_grobs
  panel_grob$children[[grep_grob(panel_grob$children, "venn")]] <- venn_grob
  mygrobs[[grep_grob(mygrobs, "panel")]] <- panel_grob
  p2$grobs <- mygrobs
  grid::grid.newpage()
  grid::grid.draw(p2)
  
}