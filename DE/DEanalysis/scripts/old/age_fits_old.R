
dir.create("plots/ctl_age")
dir.create("plots/fnf_age")

for (clust in 1:2){
  
  if (clust == 1){
    median_color <- "#267C71"
  } else if (clust == 2){
    median_color <- "#CC5E0C"
  } else if (clust == 3){
    median_color  <- "#E0AB07"
  }
  
  ctl_cluster <- all_age_ctl_lrt_data %>%
    filter(cluster == clust)
  
  # Make individual plots for each gene
  for (gene in unique(ctl_cluster$gene_id)){
    ctl_cluster_gene <- ctl_cluster %>%
      filter(gene_id == gene)
    geneSymbol <- unique(ctl_cluster_gene$symbol)
    ggplot(ctl_cluster_gene, aes(x = Age, y = logcount)) +
      geom_point(color = "grey50") +
      geom_line(aes(y = logmu), lwd = 0.75, color = median_color) +
      scale_y_continuous(name = "Log Expression") +
      theme_custom_scatterplot() +
      ggtitle(geneSymbol)
    ggsave(filename = paste0("plots/ctl_age/ctl_cluster", clust, "_", geneSymbol, ".pdf"),
           width = 7, height = 8, units = "in")
  }
  
  # Make entire count cluster plot
  ggplot(ctl_cluster, aes(x = Age, y = logcount, color = symbol)) +
    geom_line(lwd = 0.25) +
    scale_color_manual(values = rep("grey75", nrow(ctl_cluster))) +
    stat_summary(geom = "line", fun = "mean", color = median_color, lwd = 0.75) +
    scale_y_continuous(name = "Log Expression Counts") +
    theme_custom_scatterplot() +
    theme(legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5),
          legend.position = "none",
          axis.text = element_text(size = 6, color = "grey50"),
          axis.line = element_line(color = "grey50"))
  #assign(paste0("ctl_age_cluster", clust), ctl_plot, envir = globalenv())
  ggsave(filename = paste0("plots/ctl_age_cluster", clust, "_counts.pdf"),
         width = 7, height = 8, units = "in")
  
  # Make entire splinefit cluster plot
  ggplot(ctl_cluster, aes(x = Age, y = logmu, color = symbol)) +
    geom_line(lwd = 0.25) +
    scale_color_manual(values = rep("grey75", nrow(ctl_cluster))) +
    stat_summary(geom = "line", fun = "mean", color = median_color, lwd = 0.75) +
    scale_y_continuous(name = "Log Expression Fit") +
    theme_custom_scatterplot() +
    theme(legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5),
          legend.position = "none",
          axis.text = element_text(size = 6, color = "grey50"),
          axis.line = element_line(color = "grey50"))
  
  ggsave(filename = paste0("plots/ctl_age_cluster", clust, "_splines.pdf"),
         width = 7, height = 8, units = "in")
  
  fnf_cluster <- all_age_fnf_lrt_data %>%
    filter(cluster == clust)
  
  # Make individual plots for each gene
  for (gene in unique(fnf_cluster$gene_id)){
    fnf_cluster_gene <- fnf_cluster %>%
      filter(gene_id == gene)
    geneSymbol <- unique(fnf_cluster_gene$symbol)
    ggplot(fnf_cluster_gene, aes(x = Age, y = logcount)) +
      geom_point(color = "grey50") +
      geom_line(aes(y = logmu), lwd = 0.75, color = median_color) +
      scale_y_continuous(name = "Log Expression") +
      theme_custom_scatterplot() +
      ggtitle(geneSymbol)
    ggsave(filename = paste0("plots/fnf_age/fnf_cluster", clust, "_", geneSymbol, ".pdf"),
           width = 7, height = 8, units = "in")
  }
  
  
  # Make entire count cluster plot
  ggplot(fnf_cluster, aes(x = Age, y = logcount, color = symbol)) +
    geom_line(lwd = 0.25) +
    scale_color_manual(values = rep("grey75", nrow(ctl_cluster))) +
    stat_summary(geom = "line", fun = "mean", color = median_color, lwd = 0.75) +
    scale_y_continuous(name = "Log Expression Counts") +
    theme_custom_scatterplot() +
    theme(legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5),
          legend.position = "none",
          axis.text = element_text(size = 6, color = "grey50"),
          axis.line = element_line(color = "grey50"))
  
  ggsave(filename = paste0("plots/fnf_age_cluster", clust, "_counts.pdf"),
         width = 7, height = 8, units = "in")
  
  
  # Make entire splinefit cluster plot
  ggplot(fnf_cluster, aes(x = Age, y = logmu, color = symbol)) +
    geom_line(lwd = 0.25) +
    scale_color_manual(values = rep("grey75", nrow(ctl_cluster))) +
    stat_summary(geom = "line", fun = "mean", color = median_color, lwd = 0.75) +
    scale_y_continuous(name = "Log Expression Fit") +
    theme_custom_scatterplot() +
    theme(legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5),
          legend.position = "none",
          axis.text = element_text(size = 6, color = "grey50"),
          axis.line = element_line(color = "grey50"))
  
  ggsave(filename = paste0("plots/fnf_age_cluster", clust, "_splines.pdf"),
         width = 7, height = 8, units = "in")
  
  # assign(paste0("fnf_age_cluster", clust), fnf_plot, envir = globalenv())
  # ggsave(filename = paste0("plots/fnf_age_cluster", clust, ".pdf"),
  #        width = 7, height = 8, units = "in")
  
}



# Venn diagram of age cluster genes ---------------------------------------
ctl_genes <- unique(all_age_ctl_lrt_data$symbol)
fnf_genes <- unique(all_age_fnf_lrt_data$symbol)


fnf_onlygenes <- fnf_genes[which(!fnf_genes %in% ctl_genes)]
fnf_onlygenes <- all_age_fnf_lrt_data %>%
  filter(symbol %in% fnf_onlygenes)


ggplot(fnf_onlygenes, aes(x = Age, y = logmu, color = symbol)) +
  geom_line()

