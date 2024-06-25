library(ggplot2)
library(DESeq2)
library(tidyverse)
library(grid)
library(colorspace)
library(RColorBrewer)
library(ggtext)
library(extrafont)
library(showtext)

# # Add fonts ---------------------------------------------------------------
# 
# font_add(family = "Aktiv Grotesk", regular = "/work/users/n/e/nekramer/Aktiv Grotesk/OTF/AktivGrotesk-Regular.otf",
#          bold = "/work/users/n/e/nekramer/Aktiv Grotesk/OTF/AktivGrotesk-Bold.otf",
#         italic = "/work/users/n/e/nekramer/Aktiv Grotesk/OTF/AktivGrotesk-Italic.otf",
#         symbol =  "/work/users/n/e/nekramer/Aktiv.Grotesk/Aktiv Grotesk/OTF/AktivGrotesk-XBold.otf")
# # 
# font_add(family = "Aktiv Grotesk X", regular = "/work/users/n/e/nekramer/Aktiv Grotesk/OTF/AktivGrotesk-Regular.otf" ,
#          bold = "/work/users/n/e/nekramer/Aktiv.Grotesk/Aktiv Grotesk/OTF/AktivGrotesk-Hairline.otf")


# Colors ------------------------------------------------------------------

conditionColors <- c("CTL" = "#B8B8B8", "FNF" = "#4A4A4A")
ageColors <- sequential_hcl(n = 6, palette = "Mint")
ancestryColors <- c("#C74A53", "#EFB06E", "#5BAD58", "#177F97", "#704D9E")
sexColors <- c("F" = "#DD8492", "M" = "#4788BA")

heatmapColors <- colorRampPalette(c("#4FAFE1", "black", "#F5D348"))(7)
log2fcColors <- c("+" = "#FBBE67", "-" = "#78A1Cd")

ageClusterColors <- c("-" = "#009E8E", "+" = "#AC7F21")

yl_gn_bu <- brewer.pal(n = 9, name = "YlGnBu")

# Themes ---------------------------------------------------------------

theme_custom_general <- function(){
  
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(), 
          axis.ticks.x = element_blank(),
          panel.background = element_rect(color = "transparent", fill = "transparent"),
          plot.background = element_rect(color = "transparent", fill = "transparent"),
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
        panel.background = element_rect(color = "transparent", fill = "transparent"),
        plot.background = element_rect(color = "transparent", fill = "transparent"),
        axis.line = element_line(color = "black", 
                                   linewidth = 0.25),
        axis.ticks.length = unit(-0.1, "cm"),
        axis.ticks = element_line(color = "black",
                                    linewidth = 0.25),
        text = element_text(family = "Helvetica"),
        legend.key = element_blank(),
        legend.background = element_blank())
  
}


# Venn diagram fonts ------------------------------------------------------

# Function to change the fonts in a ggplot venn diagram
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
  return(p2)
  
}






# Horizontal eQTL boxplots ------------------------------------------------

create_eqtl_boxplot_horizontal <- function(data, group = NA, 
                                           stat_loc = "bottom_right",
                                           rsID_loc = "top",
                                           highlightAllele = NA,
                                           condition_labs = "top_left"){
  
  # Filter data for data from that grouping, if applicable,
  # and set genotype factors based on gt_order
  if (!is.na(group)){
    boxplot_data <- data |> 
      filter(grouping == group) |> 
      arrange(gt_order) |> 
      mutate(gt_GT_alleles = fct_inorder(gt_GT_alleles))
  } else {
    boxplot_data <- data |> 
      arrange(gt_order) |> 
      mutate(gt_GT_alleles = fct_inorder(gt_GT_alleles))
  }
  
  # Reformat gt_GT_alleles text with bolding based on highlightAllele param
  if (!is.na(highlightAllele)){
    alleles <- unlist(str_split(boxplot_data |> pull(variantID) |> unique(), ":"))[3:4]
    minor_allele <- boxplot_data |> pull(MA) |> unique()
    major_allele <- alleles[which(alleles != minor_allele)]
    # Put * for bolding around allele and reset gt factor order
    if (highlightAllele == "minor"){
      boxplot_data <- boxplot_data |> 
        mutate(gt_GT_alleles = gsub(minor_allele,
                                    paste0("**", minor_allele, "**"),
                                    gt_GT_alleles)) |> 
        arrange(gt_order) |> 
        mutate(gt_GT_alleles = fct_inorder(gt_GT_alleles))
    } else if (highlightAllele == "major"){
      
      boxplot_data <- boxplot_data |> 
        mutate(gt_GT_alleles = gsub(major_allele,
                                    paste0("**", major_allele, "**"),
                                    gt_GT_alleles)) |> 
        arrange(gt_order) |> 
        mutate(gt_GT_alleles = fct_inorder(gt_GT_alleles))
    } else if (highlightAllele == "risk") {
      risk_allele <- boxplot_data |> pull(risk_allele) |> unique()
      
      if (risk_allele == minor_allele){
        boxplot_data <- boxplot_data |> 
          mutate(gt_GT_alleles = gsub(minor_allele,
                                      paste0("**", minor_allele, "**"),
                                      gt_GT_alleles)) |> 
          arrange(gt_order) |> 
          mutate(gt_GT_alleles = fct_inorder(gt_GT_alleles))
      } else if (risk_allele == major_allele){
        boxplot_data <- boxplot_data |> 
          mutate(gt_GT_alleles = gsub(major_allele,
                                      paste0("**", major_allele, "**"),
                                      gt_GT_alleles)) |> 
          arrange(gt_order) |> 
          mutate(gt_GT_alleles = fct_inorder(gt_GT_alleles))
      }
      
    }
  }
  
  # Grab rsID and eGene name for labeling
  rsid <- unique(boxplot_data$rsID)
  eGene_name <- unique(boxplot_data$gene_symbol)
  
  # Create dummy data to label betas and pvals 
  if (stat_loc == "bottom_right"){
    stat_gt_GT_alleles <- 3.25
    stat_beta_expression <- -3.35
    stat_nom_expression <- -3.5
    ylim <- c(-3.5, 3.5)
    stat_hjust <- 1
  } else if (stat_loc == "top_right"){
    stat_gt_GT_alleles  <- 3.25
    stat_beta_expression <- 3.5
    stat_nom_expression <- 3.35
    ylim <- c(-3, 3.5)
    stat_hjust <- 1
  } else if (stat_loc == "top") {
    stat_gt_GT_alleles <- 2
    stat_beta_expression <- 3.5
    stat_nom_expression <- 3.35
    ylim <- c(-3, 3.5)
    stat_hjust <- 0.5
  }
  
  stat_data <- unique(boxplot_data[,c("beta_pbs", "beta_fnf", "nompval_pbs", "nompval_fnf")]) |> 
    pivot_longer(everything()) |> 
    separate_wider_delim(cols = "name", delim = "_", names = c("stat", "Condition")) |> 
    mutate(Condition = toupper(Condition)) |> 
    mutate(Condition = ifelse(Condition == "FNF", "FN-f", Condition)) |> 
    pivot_wider(names_from = stat, values_from = value) |> 
    mutate(Condition = factor(Condition, levels = c("PBS", "FN-f"))) |> 
    mutate(gt_GT_alleles = stat_gt_GT_alleles,
           beta_expression = stat_beta_expression,
           nom_expression = stat_nom_expression)
  
  
  
  # Create main boxplot
  eqtl_boxplot <- ggplot(boxplot_data, aes(x = gt_GT_alleles, y = expression)) +
    geom_hline(yintercept = 0, lty = 2, linewidth = 0.2) +
    geom_boxplot(aes(color = Condition, fill = Condition), outlier.shape = NA) +
    geom_jitter(position = position_jitter(width = 0.1), color = "grey25", size = 0.25) +
    # beta labels
    geom_text(data = stat_data, aes(y = beta_expression, label = paste0("Beta = ", signif(beta, digits = 3))),
              family = "Helvetica", vjust = 0, size = 1.75, hjust = stat_hjust) +
    # p-value labels
    geom_text(data = stat_data, aes(y = nom_expression, label = paste0("p-value = ", signif(nompval, digits = 3))),
              family = "Helvetica", vjust = 1, size = 1.75, hjust = stat_hjust) +
    
    scale_color_manual(values = c("#48617b", "#97723e")) +
    scale_fill_manual(values = c("#78A1Cd", "#FBBE67"))  + 
    facet_wrap(vars(Condition), nrow = 1, strip.position = "left") +
    scale_y_continuous(name = paste0("**", eGene_name, "**", " normalized expression"),
                       limits = ylim, breaks = seq(-3, 3, 1)) +
    theme(legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(), 
          axis.ticks.x = element_blank(),
          axis.line.y = element_line(color = "black", 
                                     linewidth = 0.2),
          axis.ticks.length.y = unit(-0.1, "cm"),
          axis.ticks.y = element_line(color = "black",
                                      linewidth = 0.2),
          panel.background = element_rect(fill = 'transparent', color = "transparent"),
          plot.background = element_rect(fill = 'transparent', color = "transparent"),
          strip.background = element_blank(),
          strip.placement = "outside",
          strip.text.y.left = element_blank(),
          strip.text = element_text(color = "black", size = 11),
          axis.text = element_text(color = "black"),
          axis.text.x = element_markdown(margin = margin(t=-8)),
          text = element_text(family = "Helvetica", size = 8),
          axis.title.y = element_markdown(size = 8, margin = margin(r = -2)),
          plot.margin = margin(l = 0, r = 0, t = 10, b = 10),
          panel.spacing = unit(0, "cm")) +
    coord_cartesian(clip = "off")
  
  # Add rsID label to boxplot object
  if (rsID_loc == "top"){
    eqtl_boxplot <- eqtl_boxplot +
      ggtitle(rsid) +
      theme(axis.title.x = element_blank(),
            plot.title = element_text(hjust = 0.5, size = 7,
                                      margin = margin(b = -1)))
  } else if (rsID_loc == "bottom"){
    eqtl_boxplot <- eqtl_boxplot +
      xlab(rsid) +
      theme(axis.title.x = element_text(size = 7, vjust = 6),
            plot.title = element_text(hjust = 0.5, size = 10, face = "bold")) +
      ggtitle(group)
  }
  
  # Condition labels
  if (!is.na(condition_labs)){
    if (condition_labs == "top_left") {
      
      # Create dummy data to add PBS and FN-f labels to top
      condition_data <- tibble(Condition = c("PBS", "FN-f"),
                               gt_GT_alleles = 0.5,
                               expression = 3.4)
      condition_data$Condition <- factor(condition_data$Condition, levels = c("PBS", "FN-f"))
      eqtl_boxplot <- eqtl_boxplot +
        # Condition labels
        geom_text(data = condition_data, aes(y = expression, label = Condition, color = Condition), 
                  family = "Helvetica", size = 3, hjust = 0)
      
      
    } else if (condition_labs == "bottom_left") {
      # Create dummy data to add PBS and FN-f labels to bottom
      condition_data <- tibble(Condition = c("PBS", "FN-f"),
                               gt_GT_alleles = 0.5,
                               expression = -3)
      condition_data$Condition <- factor(condition_data$Condition, levels = c("PBS", "FN-f"))
      eqtl_boxplot <- eqtl_boxplot +
        # Condition labels
        geom_text(data = condition_data, aes(y = expression, label = Condition, color = Condition), 
                  family = "Helvetica", size = 2.5, hjust = 0)
      
    } else if (condition_labs == "bottom") {
      condition_data <- tibble(Condition = c("PBS", "FN-f"),
                               gt_GT_alleles = 2,
                               expression = -3)
      condition_data$Condition <- factor(condition_data$Condition, levels = c("PBS", "FN-f"))
      eqtl_boxplot <- eqtl_boxplot +
        # Condition labels
        geom_text(data = condition_data, aes(y = expression, label = Condition, color = Condition), 
                  family = "Helvetica", size = 2.5, hjust = 0.5, vjust = 0)
    }
  }
  
  

  return(eqtl_boxplot)
  
}
