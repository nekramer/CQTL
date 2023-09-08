plotGOKEGG_separate <- function(data, cat, type){
  
  if (cat == "up"){
    label <- "Upregulated Genes"
    axiscolor <- "#FBBE67"
  } else {
    label <- "Downregulated Genes"
    axiscolor <- "#78A1Cd"
  }
  
  if (type == "GO"){
    yName <- "parentTerm"
  } else {
    yName <- "Term"
  }
  
  GOKEGGplot <- ggplot(data, aes(x = category, y = get(yName))) +
    geom_tile(aes(fill =  `-log10pval`)) +
    scale_y_discrete(position = "right") +
    scale_x_discrete(labels = c(label)) +
    scale_fill_gradientn(name = "-log10 p-value", 
                         
                         colors = brewer.pal(n = 7, name = "YlOrRd"),
                         limits = c(0, 
                                    ceiling(max(data$`-log10pval`))),
                         breaks = c(0,
                                    ceiling(max(data$`-log10pval`)))) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(), 
          axis.ticks.x = element_blank(),
          panel.background = element_blank(), 
          axis.line.y = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.text.y = element_text(color = "black", size = 12),
          text = element_text(family = "Helvetica"),
          legend.position = "top",
          legend.justification = 0.5,
          legend.title = element_text(size = 10, vjust = -5.5),
          legend.text = element_text(size = 10, vjust = 0),
          axis.text.x = element_text(color = axiscolor, face = "bold", size = 12,
                                     hjust = 0.5)) +
    guides(fill = guide_colorbar(title.position = "top",
                                 title.hjust = 0.5,
                                 ticks = FALSE,
                                 label.position = "top",
                                 barwidth = unit(1.75, "in"))) +
    coord_fixed(ratio = 0.5)
  
  ggsave(filename = paste0("plots/", cat, "sig_", type, ".pdf"),
         plot = GOKEGGplot,
         width = 8, height = 10, units = "in",
         bg = "transparent")
  
  save(GOKEGGplot, file = paste0("data/", cat, "sig_", type, ".rda"))
  
}

plotGOKEGG_overlap <- function(data, type){
  
  if (type == "GO"){
    yName <- "term"
  } else {
    yName <- "Term"
  }
  
  upregsubset <- ggplot(data, mapping = aes(x = category, y = get(yName))) +
    geom_tile(aes(fill = `-log10pval`)) +
    scale_y_discrete(position = "right") +
    scale_x_discrete(labels = c(str_wrap("Upregulated Genes", width = 5),
                                str_wrap("Downregulated Genes", width = 5))) +
    scale_fill_gradientn(name = "-log10 p-value", 
                         
                         colors = brewer.pal(n = 7, name = "YlOrRd"),
                         limits = c(0, 
                                    ceiling(max(data$`-log10pval`))),
                         breaks = c(0,
                                    ceiling(max(data$`-log10pval`)))) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(), 
          axis.ticks.x = element_blank(),
          panel.background = element_blank(), 
          axis.line.y = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.text.y = element_text(color = "black", size = 12),
          text = element_text(family = "Helvetica"),
          legend.position = "top",
          legend.justification = 0.5,
          legend.title = element_text(size = 10, vjust = -5.5),
          legend.text = element_text(size = 10, vjust = 0),
          axis.text.x = element_text(color = c("#FBBE67", "#78A1Cd"), face = "bold", size = 10,
                                     hjust = 0.5)) +
    guides(fill = guide_colorbar(title.position = "top",
                                 title.hjust = 0.5,
                                 ticks = FALSE,
                                 label.position = "top",
                                 barwidth = unit(1.75, "in"))) +
    coord_fixed(ratio = 0.25)
  
  ggsave(filename = paste0("plots/", type, "_upregsubset.pdf"), width = 8, height = 10, 
         units = "in")
  
  save(upregsubset, file = paste0("data/", type, "_upregsubset.rda"))
  
}

plotGOKEGG_stacked <- function(upData, downData, type, scale = "0"){
  
  if (type == "GO"){
    title <- "GO Terms"
    yName <- "parentTerm"
  } else {
    title <- "KEGG Pathways"
    yName <- "Term"
  }
  
  if (scale == "0"){
    minlim_up <- 0
    minlim_down <- 0
  } else if (scale == "min"){
    minlim_up <- floor(min(upData$`-log10pval`))
    minlim_down <- floor(min(downData$`-log10pval`))
  }
  
  
  up_heatmap <- ggplot(data = upData, aes(x = category, y = get(yName))) +
    geom_tile(aes(fill = `-log10pval`)) +
    scale_fill_gradientn(name = "-log10 p-value", 
                         
                         colors = c("white", "#FBBE67"),
                         limits = c(minlim_up, 
                                    ceiling(max(upData$`-log10pval`))),
                         breaks = c(minlim_up,
                                    ceiling(max(upData$`-log10pval`)))) +
    scale_y_discrete(position = "right") +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(), 
          axis.ticks.x = element_blank(),
          panel.background = element_blank(), 
          axis.line.y = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.text.y = element_text(color = "black", size = 12),
          text = element_text(family = "Helvetica"),
          axis.text.x = element_blank(),
          plot.title = element_text(hjust = 0.75, 
                                    vjust = -2,
                                    face = "bold", 
                                    color = "#FBBE67",
                                    size = 11),
          plot.title.position = "plot",
          legend.position = "left",
          plot.margin = margin(b = 0),
          legend.title = element_text(angle = 90, size = 10),
          legend.text = element_text(margin = margin(l = -19),
                                     hjust = 0.25,
                                     vjust = c(1.5, -0.5))) +
    guides(fill = guide_colorbar(title.position = "left",
                                 title.hjust = 0.5,
                                 ticks = FALSE)) +
    ggtitle("Upregulated Genes") +
    coord_fixed(ratio = 0.75)
  
  
  
  down_heatmap <- ggplot(data = downData, aes(x = category, y = get(yName))) +
    geom_tile(aes(fill = `-log10pval`)) +
    scale_fill_gradientn(name = "-log10 p-value", 
                         
                         colors = c("white", "#78A1Cd"),
                         limits = c(minlim_down, 
                                    ceiling(max(downData$`-log10pval`))),
                         breaks = c(minlim_down,
                                    ceiling(max(downData$`-log10pval`)))) +
    scale_y_discrete(position = "right") +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(), 
          axis.ticks.x = element_blank(),
          panel.background = element_blank(), 
          axis.line.y = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.text.y = element_text(color = "black", size = 12),
          text = element_text(family = "Helvetica"),
          axis.text.x = element_blank(),
          plot.title = element_text(hjust = 0.75,
                                    vjust = -2,
                                    face = "bold", 
                                    color = "#78A1Cd",
                                    size = 11),
          plot.title.position = "plot",
          legend.position = "left",
          plot.margin = margin(t = 0),
          legend.title = element_text(angle = 90, size = 10),
          legend.text = element_text(margin = margin(l = -19),
                                     hjust = 0.25,
                                     vjust = c(1.5, -0.5))) +
    guides(fill = guide_colorbar(title.position = "left",
                                 title.hjust = 0.5,
                                 ticks = FALSE,
                                 reverse = TRUE)) +
    ggtitle("Downregulated Genes") +
    coord_fixed(ratio = 0.75)
  
  
  stacked_heatmap <- up_heatmap + down_heatmap + 
    plot_layout(ncol = 1) +
    plot_annotation(title = title, 
                    theme = theme(text = element_text(family = "Helvetica",
                                                      face = "bold")))
  
  
  ggsave(filename = paste0("plots/", type, "stacked_heatmap_", scale,".pdf"),
         plot = stacked_heatmap,
         width = 7, height = 9, units = "in",
         bg = "transparent")
  
  #save(stacked_heatmap, file = paste0("data/", type, "stacked_heatmap.rda"))
}