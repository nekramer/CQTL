## Vertical genes
for (motif in unique(upsig_knownmotifs$Name)){
  motif_genes <- upsig_knownmotifs %>%
    filter(Name == motif)
  
  motif_gene_expression <- list()
  
  for (gene in motif_genes$Gene){
    # Grab ENSEMBL ID
    geneSymbol <- mapIds(org.Hs.eg.db, keys = gene, keytype = "SYMBOL",
                         column = "ENSEMBL")
    
    # Get counts
    countData <- plotCounts(dds, gene = geneSymbol, 
                            intgroup = "Condition", 
                            returnData = TRUE) %>%
      mutate(Gene = gene)
    
    
    motif_gene_expression[[gene]] <- countData
  }
  
  motif_gene_expression <- bind_rows(motif_gene_expression) %>%
    remove_rownames()
  
  # Change text of NFkb to include symbol with unicode
  if (motif == "NFkB-p65-Rel"){
    motif <- "NF-\u03BAB"
  }
  
  motifP <- unique(motif_genes$`P-value`)
  # Create inset motif image with title and p-value and make grob 
  motifImg <- pictureGrob(readPicture(rawToChar(rsvg::rsvg_svg(unique(motif_genes %>% pull(motifLogo))))))
  motifGrab <- grid.grabExpr(expr = {
    grid.newpage()
    grid.draw(motifImg)
    grid.text(label = motif, x = 0.5, y = 0.35, 
              gp = gpar(fontfamily = "Helvetica",
                        fontsize = 14))
    grid.text(label = paste0("p-value: ", motifP), 
              x = 0.5, y = 0.25, 
              gp = gpar(fontfamily = "Helvetica",
                        fontsize = 12))
  })
  
  geneLimits <- unlist(str_split(paste(unique(motif_gene_expression$Gene), "skip"), " "))
  geneLimits <- geneLimits[-length(geneLimits)]
  boxplots <- ggplot(motif_gene_expression, aes(x = count, y = Gene, fill = Condition)) +
    geom_boxplot() +
    scale_x_continuous(name = "Normalized counts", expand = c(0, 0),
                       #limits = c(0, round(max(motif_gene_expression$count)) + 1000),
                       trans = "log10") +
    scale_y_discrete(breaks = unique(motif_gene_expression$Gene),
                     limits = geneLimits) +
    scale_fill_manual(values = directionColors) +
    theme(legend.position = "none",
          axis.text.x = element_text(size = 12),
          axis.text.y = element_markdown(family = "Helvetica", size = 12),
          axis.title = element_text(size = 14),
          axis.title.y = element_blank()) +
    theme_custom_scatterplot()
  # grDevices::cairo_pdf(paste0("plots/upMotif_", motif, "_horizontal.pdf"), 
  #                      width = 5, height = 5)
  wrap_plots(A = motifGrab, B = boxplots, widths = c(0.5, 1), design = "AB
                                                                        #B")
  
  # dev.off()
  ggsave(paste0("plots/upMotif_", motif, "_horizontal.pdf"),width = 5, height = 5,
         device = cairo_pdf)
  
}


## Horizontal genes
for (motif in unique(upsig_knownmotifs$Name)){
  
  motif_genes <- upsig_knownmotifs %>%
    filter(Name == motif)
  
  motif_gene_expression <- list()
  
  for (gene in motif_genes$Gene){
    # Grab ENSEMBL ID
    geneSymbol <- mapIds(org.Hs.eg.db, keys = gene, keytype = "SYMBOL",
                         column = "ENSEMBL")
    
    # Get counts
    countData <- plotCounts(dds, gene = geneSymbol, 
                            intgroup = "Condition", 
                            returnData = TRUE) %>%
      mutate(Gene = gene) %>%
      mutate(sig = ifelse(Gene %in% sig_upgenes$symbol, "*", ""))
    
    
    motif_gene_expression[[gene]] <- countData
  }
  
  motif_gene_expression <- bind_rows(motif_gene_expression) %>%
    remove_rownames()  %>%
    mutate(geneSig = ifelse(sig == "*", paste0(Gene,"^", sig, "^"), Gene))
  
  # Relevel FOSL2 so high expression doesn't block motif image
  if (motif == "Fos"){
    motif_gene_expression$geneSig <- factor(motif_gene_expression$geneSig)
    motif_gene_expression$geneSig <-relevel(motif_gene_expression$geneSig, "FOSL2")
  }
  
  # Change text of NFkb to include symbol with unicode
  if (motif == "NFkB-p65-Rel"){
    motif <- "NF-\u03BAB"
  }
  
  motifP <- unique(motif_genes$`P-value`)
  # Create inset motif image with title and p-value and make grob 
  motifImg <- pictureGrob(readPicture(rawToChar(rsvg::rsvg_svg(unique(motif_genes %>% pull(motifLogo))))))
  motifGrab <- grid.grabExpr(expr = {
    grid.newpage()
    grid.draw(motifImg)
    grid.text(label = motif, x = 0.95, y = 0.25, just = c("right", "bottom"), 
              gp = gpar(fontfamily = "Helvetica",
                        fontsize = 14))
    grid.text(label = paste0("p-value: ", motifP), 
              x = 0.95, y = 0.15, just = c("right", "bottom"), 
              gp = gpar(fontfamily = "Helvetica",
                        fontsize = 12))
  })
  
  # Use cairo graphic device for appropriate rendering of "kappa" unicode character
  grDevices::cairo_pdf(paste0("plots/upMotif_", motif, ".pdf"), 
                       width = 6, height = 6)
  # Plot boxplots for genes based on motifs
  ggplot(motif_gene_expression, aes(x = geneSig, y = count, fill = Condition)) +
    geom_boxplot() +
    scale_x_discrete(name = "TF Gene") +
    scale_y_continuous(name = "Normalized counts", expand = c(0, 0),
                       limits = c(0, round(max(motif_gene_expression$count)) + 1000)) +
    scale_fill_manual(values = directionColors) +
    theme(legend.position = "none",
          axis.text.y = element_text(size = 12),
          axis.text.x = element_markdown(family = "Helvetica", size = 12),
          axis.title = element_text(size = 14)) +
    theme_custom_scatterplot() +
    inset_element(p = motifGrab, left = 0.7, bottom = 0.8, right = 1, top = 1.1,
                  align_to = "full", on_top = FALSE)
  
  dev.off()
  
}
