library(ComplexHeatmap)
library(DESeq2)
library(tidyverse)
library(circlize)
library(RColorBrewer)
library(plotgardener)
library(grid)
library(rrvgo)
library(ggnewscale)
library(patchwork)
library(rsvg)
library(grImport2)
library(ggrepel)
library(org.Hs.eg.db)
source("../plotting_utils.R")
source("../utils.R")
library(gridtext)

# Functions ---------------------------------------------------------------

# Function to make grid-based boxplots of log2FoldChanges with highlighted genes
oa_boxplot_plot <- function(data, up_wilcox_test, down_wilcox_test, x, y, width, height,
                            default.units = "inches", just = c("left", "top")) {
  
  ## Create a function to draw boxplot with distribution
  makeBoxplot <- function(vec,
                          fill,
                          col,
                          yscale,
                          xpos = 0.2, 
                          width = 0.25, 
                          alpha = 0.6){
    
    ## Calculate median, q1, q3, 1.5*IQR
    boxplot_median <- median(vec)
    boxplot_q1 <- quantile(vec, 0.25)
    boxplot_q3 <- quantile(vec, 0.75)
    boxplot_IQR <- boxplot_q3 - boxplot_q1
    
    ## Define viewport
    
    vp <- viewport(x = unit(xpos, "npc"), 
                   y = unit(0.5, "npc"),
                   width = width, 
                   height = 1, just = c("left", "center"),
                   yscale = yscale,
                   clip = "off",
                   name = "boxplot")
    
    # Create boxplot background box
    
    boxplot_box <- polygonGrob(
      # bottom left, top left, top right, bottom right
      x = c(unit(0, "npc"), unit(0, "npc"), unit(1, "npc"), unit(1, "npc")),
      y = c(unit(boxplot_q1, "native"), unit(boxplot_q3, "native"), 
            unit(boxplot_q3, "native"), unit(boxplot_q1, "native")),
      default.units = "native",
      gp = gpar(fill = alpha(colour = fill, alpha = alpha),
                col = col, lwd = 1)
    )
    
    ## Boxplot distribution lines
    # Median
    median_line <- segmentsGrob(x0 = 0,
                                y0 = unit(boxplot_median, "native"),
                                x1 = 1, y1 = unit(boxplot_median, "native"),
                                gp = gpar(lwd = 3, lineend = "butt", col = col))
    
    # Q1
    q1_line <- segmentsGrob(x0 = 0,
                            y0 = unit(boxplot_q1, "native"),
                            x1 = 1, y1 = unit(boxplot_q1, "native"),
                            gp = gpar(lwd = 1, lineend = "square", col = col))
    
    # Q3
    q3_line <- segmentsGrob(x0 = 0,
                            y0 = unit(boxplot_q3, "native"),
                            x1 = 1, y1 = unit(boxplot_q3, "native"),
                            gp = gpar(lwd = 1, lineend = "square", col = col))
    
    
    # Q3 + 1.5*IQR
    q3_iqr_line <- segmentsGrob(x0 = 0.25,
                                y0 = unit(boxplot_q3 + 1.5*boxplot_IQR, "native"),
                                x1 = 0.75,
                                y1 = unit(boxplot_q3 + 1.5*boxplot_IQR, "native"),
                                gp = gpar(lwd = 2, lineend = "square", col = col))
    
    # Q1 - 1.5*IQR
    q1_iqr_line <- segmentsGrob(x0 = 0.25,
                                y0 = unit(boxplot_q1 - 1.5*boxplot_IQR, "native"),
                                x1 = 0.75,
                                y1 = unit(boxplot_q1 - 1.5*boxplot_IQR, "native"),
                                gp = gpar(lwd = 2, lineend = "square", col = col))
    
    # Vertical line
    v_line <- segmentsGrob(x0 = 0.5, y0 = unit(boxplot_q1 - 1.5*boxplot_IQR, "native"),
                           x1 = 0.5, y1 = unit(boxplot_q3 + 1.5*boxplot_IQR, "native"),
                           gp = gpar(lwd = 1, lineend = "square", col = col))
    
    boxplot_gtree <- gTree(vp = vp, children = gList(boxplot_box, 
                                                     median_line,
                                                     q1_line,
                                                     q3_line,
                                                     q3_iqr_line,
                                                     q1_iqr_line,
                                                     v_line))
    return(boxplot_gtree)
  }
  
  ## Create a function to add jittered points
  makePoints <- function(vec, yscale,
                         xpos = 0.2, 
                         xwidth = 0.05,
                         pch = 21, size = 0.5, 
                         fill = "black", alpha = 1, col = NA, 
                         jitter = TRUE, jwidth = 0.5){
    get_xval_diff <- function(val){
      if (as.numeric(val) < 0) {
        xval_diff <- 1 - abs(as.numeric(val))
      } else {
        xval_diff <- 1 + as.numeric(val)
      }
    }
    
    ## Toggle jittered xposition
    if (jitter){
      set.seed(123)
      xvals <- unit(runif(length(vec), -jwidth, jwidth), "npc")
    } else {
      xvals <- unit(rep(0, length(vec)), "npc")
    }
    
    ## Define viewport
    vp <- viewport(x = unit(xpos, "npc"), y = unit(0.5, "npc"), width = xwidth, height = 1,
                   just = c("left", "center"), yscale = yscale,
                   name = "points")
    ## Points
    points <- pointsGrob(x = xvals, y = unit(vec, "native"),
                         pch = pch, size = unit(size, "npc"), 
                         gp = gpar(fill = fill, col = col, alpha = alpha))
    
    # Create gTree
    points_gtree <- gTree(vp = vp, children = gList(points))
    
    ## Return xvalues for labels
    xvals_diff <- unlist(lapply(xvals, get_xval_diff))
    names(xvals_diff) <- names(vec)
    
    return(list(points_gtree, xvals_diff))
  }
  
  yscale <- c(-4, 6)
  
  plotVP <- viewport(x = x, y = y, width = width, height = height,
                     default.units = default.units, just = just,
                     yscale = yscale, name = "plot_area")
  
  oa_boxplot_gtree <- gTree(name = "OA_boxplots", 
                            vp = plotVP)
  
  ## Y- axis 
  ylabs <- seq(-4, 6, 2)
  
  # Axis line
  oa_boxplot_gtree <- addGrob(oa_boxplot_gtree, 
                              child = segmentsGrob(x0 = 0.01, 
                                                   y0 = 0, x1 = 0.01, y1 = 1,
                                                   gp = gpar(lwd = 0.75)))
  
  # axis text y
  oa_boxplot_gtree <- addGrob(oa_boxplot_gtree, 
                              child = textGrob(label = ylabs, 
                                               x = -0.01, 
                                               y = unit(ylabs, "native"),
                                               gp = gpar(fontfamily = "Helvetica", 
                                                         fontsize = 6)))
  # axis title y 
  oa_boxplot_gtree <- addGrob(oa_boxplot_gtree, 
                              child = richtext_grob(text = 
                                                      "log~2~(fold change)<br> in response to FN-f", 
                                                    x = -0.05, y = 0.5, rot = 90,
                                                    gp = gpar(fontfamily = "Helvetica", 
                                                              fontsize = 6)))
  # 0 line
  oa_boxplot_gtree <- addGrob(oa_boxplot_gtree, 
                              child = segmentsGrob(x0 = 0.01, 
                                                   y0 = unit(0, "native"), 
                                                   x1 = 0.575, 
                                                   y1 = unit(0, "native"), 
                                                   default.units = "npc",
                                                   gp = gpar(col = "grey25", lty = 2)))
  
  ## Calculate x-positioning for groups
  n <- 10
  xpos <- seq(1,2*n, 2)/(2*n)
  xpos <- xpos + 0.05 # shift xpos
  
  ## Define barwidth
  barwidth <- 0.08
  
  
  ## Jittered points
  # Up in OA, not genes of interest
  oa_up <- makePoints(vec = data |> 
                        filter(group == "Up in OA" & is.na(highlight)) |> 
                        pull(log2FoldChange), 
                      yscale = yscale,
                      xpos = xpos[2], 
                      size = 0.25, fill = "grey80", jwidth = 0.8)
  
  ## Add gTree to plot gTree
  oa_boxplot_gtree <- addGrob(oa_boxplot_gtree, child = oa_up[[1]])
  
  # Up in OA, genes of interest
  up_vec <- data |> filter(group == "Up in OA" & highlight == "up") 
  up_vec_data <- up_vec[["log2FoldChange"]]
  names(up_vec_data) <- up_vec[["symbol"]]
  
  oa_up_highlight <- makePoints(vec = up_vec_data, 
                                yscale = yscale,
                                xpos = xpos[2], 
                                size = 0.25, fill = darken(log2fcColors[["+"]], 0.3), jwidth = 0.8)
  
  oa_boxplot_gtree <- addGrob(oa_boxplot_gtree, child = oa_up_highlight[[1]])
  
  # Down in OA, not genes of interest
  oa_down <- makePoints(vec = data |> 
                          filter(group == "Down in OA" & is.na(highlight)) |> 
                          pull(log2FoldChange), 
                        yscale = yscale,
                        xpos = xpos[4], 
                        size = 0.25, fill = "grey80", jwidth = 0.8)
  
  oa_boxplot_gtree <- addGrob(oa_boxplot_gtree, child = oa_down[[1]])
  
  # Down in OA, genes of interest
  down_vec <- data |> filter(group == "Down in OA" & highlight == "down") 
  down_vec_data <- down_vec[["log2FoldChange"]]
  names(down_vec_data) <- down_vec[["symbol"]]
  
  oa_down_highlight <- makePoints(vec = down_vec_data, 
                                  yscale = yscale,
                                  xpos = xpos[4], 
                                  size = 0.25, 
                                  fill = darken(log2fcColors[["-"]], 0.3), jwidth = 0.8)
  
  oa_boxplot_gtree <- addGrob(oa_boxplot_gtree, child = oa_down_highlight[[1]])
  
  ## Boxplots
  up_boxplot <- makeBoxplot(data |> 
                              filter(group == "Up in OA") |> 
                              pull(log2FoldChange),
                            yscale = yscale,
                            fill = log2fcColors[["+"]],
                            col = darken(log2fcColors[["+"]], 0.3),
                            xpos = xpos[2]-0.0375, width = barwidth)
  oa_boxplot_gtree <- addGrob(oa_boxplot_gtree, child = up_boxplot)
  down_boxplot <- makeBoxplot(data |> 
                                filter(group == "Down in OA") |> 
                                pull(log2FoldChange),
                              yscale = yscale,
                              fill = log2fcColors[["-"]],
                              col = darken(log2fcColors[["-"]], 0.3),
                              xpos = xpos[4]-0.0375, width = barwidth)
  oa_boxplot_gtree <- addGrob(oa_boxplot_gtree, child = down_boxplot)
  
  # Get genes of interest
  goi <- data |> 
    filter(!is.na(highlight)) |> 
    pull(symbol)
  
  ## Upregulated
  goi_LFC_up <- data |> 
    filter(group == "Up in OA" & symbol %in% goi) |> 
    arrange(log2FoldChange) 
  
  ypos2 <- seq(0.5, 0.8, length.out = 5)
  ## Draw segments connecting y-positions (left)
  left_seg1 <- segmentsGrob(x0 =  (xpos[2] - (0.5*barwidth)) - 0.055, y0 = ypos2,
                            x1 =  (xpos[2] - (0.5*barwidth)) - 0.04, y1 = ypos2,
                            gp = gpar(lty=3, col = darken(log2fcColors[["+"]], 0.3)))
  oa_boxplot_gtree <- addGrob(oa_boxplot_gtree, child = left_seg1)
  
  left_seg2 <- segmentsGrob(x0 = (xpos[2] - (0.5*barwidth)) - 0.04, y0 = ypos2,
                            x1 = (xpos[2] - (0.5*barwidth)) - 0.02, 
                            y1 = unit(goi_LFC_up |> pull(log2FoldChange), "native"),
                            gp = gpar(lty=3, col = darken(log2fcColors[["+"]], 0.3)))
  oa_boxplot_gtree <- addGrob(oa_boxplot_gtree, child = left_seg2)
  
  left_seg3 <- segmentsGrob(x0 = (xpos[2] - (0.5*barwidth)) - 0.02, 
                            y0 = unit(goi_LFC_up  |> arrange(match(symbol, names(oa_up_highlight[[2]]))) |> pull(log2FoldChange), "native"),
                            x1 = (xpos[2] - 0.5*barwidth) + oa_up_highlight[[2]]*barwidth*0.5, 
                            y1  = unit(goi_LFC_up  |> arrange(match(symbol, names(oa_up_highlight[[2]]))) |> pull(log2FoldChange), "native"),
                            gp = gpar(lty = 3, col = darken(log2fcColors[["+"]], 0.3), lineend = "butt"))
  oa_boxplot_gtree <- addGrob(oa_boxplot_gtree, child = left_seg3)
  
  ## Up gene labels
  upGenes <- textGrob(label = goi_LFC_up |> pull(symbol), 
                      x = (xpos[2] - (0.5*barwidth)) - 0.06, y = ypos2, just = c("right", "center"),
                      gp = gpar(col = darken(log2fcColors[["+"]], 0.3), fontsize = 6,
                                fontfamily = "Helvetica"))
  oa_boxplot_gtree <- addGrob(oa_boxplot_gtree, child = upGenes)
  
  ## Downregulated
  goi_LFC_down <- data |> 
    filter(group == "Down in OA" & symbol %in% goi) |> 
    arrange(log2FoldChange) 
  
  
  ypos2 <- seq(0.05, 0.35, length.out = 5)
  
  ## Draw segments connecting y-positions (right)
  right_seg1 <- segmentsGrob(x0 = (xpos[4] + (0.5*barwidth)) + 0.02, 
                             y0 = unit(goi_LFC_down |> 
                                         arrange(match(symbol, names(oa_down_highlight[[2]]))) |> 
                                         pull(log2FoldChange), "native"),
                             x1 = (xpos[4] - 0.5*barwidth) + oa_down_highlight[[2]]*barwidth*0.55, 
                             y1  = unit(goi_LFC_down |> arrange(match(symbol, names(oa_down_highlight[[2]]))) 
                                        |> pull(log2FoldChange), "native"),
                             gp = gpar(lty = 3, col = darken(log2fcColors[["-"]], 0.3),
                                       lineend = "butt"))
  oa_boxplot_gtree <- addGrob(oa_boxplot_gtree, child = right_seg1)
  
  right_seg2 <- segmentsGrob(x0 = (xpos[4] + (0.5*barwidth)) + 0.02, 
                             y0 = unit(goi_LFC_down |> pull(log2FoldChange), "native"),
                             x1 = (xpos[4] + (0.5*barwidth)) + 0.04, 
                             y1 = ypos2,
                             gp = gpar(lty=3, col = darken(log2fcColors[["-"]], 0.3)))
  oa_boxplot_gtree <- addGrob(oa_boxplot_gtree, child = right_seg2)
  
  right_seg3 <- segmentsGrob(x0 =  (xpos[4] + (0.5*barwidth)) + 0.04, y0 = ypos2,
                             x1 =  (xpos[4] + (0.5*barwidth)) + 0.055, y1 = ypos2,
                             gp = gpar(lty=3, col = darken(log2fcColors[["-"]], 0.3)))
  oa_boxplot_gtree <- addGrob(oa_boxplot_gtree, child = right_seg3)
  
  
  ## Down gene labels
  downGenes <- textGrob(label = goi_LFC_down |> pull(symbol), 
                        x = (xpos[4] + (0.5*barwidth)) + 0.06, y = ypos2, just = c("left", "center"),
                        gp = gpar(col = darken(log2fcColors[["-"]], 0.3), fontsize = 6,
                                  fontfamily = "Helvetica"))
  oa_boxplot_gtree <- addGrob(oa_boxplot_gtree, child = downGenes)  
  
  ## X-axis labels
  xaxis_labels <- textGrob(label = c("Up in OA", "Down in OA"),
                           x = xpos[c(2, 4)], y =-0.025, just = c("center", "top"),
                           gp = gpar(lineheight = 0.9, fontfamily = "Helvetica",
                                     fontsize = 8, 
                                     col = c(darken(log2fcColors[["+"]], 0.3), darken(log2fcColors[["-"]], 0.3))))
  oa_boxplot_gtree <- addGrob(oa_boxplot_gtree, child = xaxis_labels)
  
  ## Significance stars from Wilcox testing enrichment
  # Upregulated
  if (up_wilcox_test$p.value < 0.01){
    up_test_star <- textGrob(label = "*", x = xpos[2], y = 0.95,
                             gp = gpar(fontface = "bold", fontfamily = "Helvetica",
                                       fontsize = 14))
  } else {
    up_test_star <- textGrob(label = "ns", x = xpos[2], y = 0.95,
                             gp = gpar(fontface = "bold", fontfamily = "Helvetica"))
  }
  oa_boxplot_gtree <- addGrob(oa_boxplot_gtree, child = up_test_star)
  # Downregulated
  if (down_wilcox_test$p.value < 0.01){
    down_test_star <- textGrob(label = "*", x = xpos[4], y = 0.95,
                               gp = gpar(fontface = "bold", fontfamily = "Helvetica",
                                         fontsize = 14))
  } else {
    down_test_star <- textGrob(label = "ns", x = xpos[4], y = 0.95,
                               gp = gpar(fontface = "bold", fontfamily = "Helvetica"))
  }
  oa_boxplot_gtree <- addGrob(oa_boxplot_gtree, child = down_test_star)
  
  
  grid.draw(oa_boxplot_gtree)
  
}

# Function to calculate sample-level l2fc for a gene in a count matrix
get_sample_l2fc <- function(gene, countMatrix){
  
  # Extract row of gene from countMatrix
  gene_counts <- countMatrix[gene,]
  
  # Convert to dataframe and extract donors/conditions into separate columns
  donor_gene_counts <- data.frame(gene_counts) |> 
    rownames_to_column(var = "Sample") |> 
    separate_wider_delim("Sample", 
                         delim = "_", 
                         names = c(NA, "Donor", NA, "Condition", NA, NA)) |> 
    # Group by each donor and calculate l2FC 
    group_by(Donor) |> 
    summarize(log2FC = 
                log2(gene_counts[Condition == "FNF"]/gene_counts[Condition == "CTL"])) |> 
    ungroup() |> 
    mutate(ENSEMBL = gene)
  
  return(donor_gene_counts)
}

# HEATMAP -----------------------------------------------------------------

load("data/condition_de/differential_expression_dds.rda")
sig_degenes <- 
  read_csv("data/condition_de/sig_deGenes_pval01_l2fc2.csv") |> 
  mutate(log2FC_dir = ifelse(log2FoldChange < 0, "-", "+")) |> 
  arrange(log2FC_dir)
donorSamplesheet <- read_csv("data/donorSamplesheet.csv") |> 
  mutate(Race = replace_na(Race, "Unknown")) |> 
  dplyr::select(Donor, Sex, Age) |> 
  # Read in and join ancestries determined through genotyping pca
  left_join(read_csv("data/CQTL_COA_01_GDA8_COA2_01_COA3_01_GDA8_COA4_COA5_COA6_COA7_predictedAncestry.csv") |> 
              separate_wider_delim(cols = "Donor", delim = "_", 
                                   names = c(NA, "Donor", NA), too_many = "drop"), by = "Donor")

# Normalized counts
dds_norm <- vst(dds)

# Grab significant genes
normCounts <- assay(dds_norm)[sig_degenes$gene_id,] |>  
  as.data.frame() 
normCounts <- normCounts[match(sig_degenes$gene_id, rownames(normCounts)),]

# Sort columns into CTL and FNF
normCounts <- normCounts %>%
  # Sort columns into CTL and FNF
  dplyr::select(contains(c("CTL", "FNF")))

# Scale counts
mat_scaled <- t(apply(normCounts, 1, scale))
colnames(mat_scaled) <- colnames(normCounts)

# Condition, Age, Sex, Race Clusters
annotations <- as.data.frame(colData(dds)[,c("Condition", "Donor")]) |> 
  rownames_to_column(var = "Sample") |> 
  # Put in same order as matrix
  arrange(match(Sample, colnames(normCounts))) |> 
  dplyr::select(-Sample)
colnames(annotations) <- c("Condition", "Donor")
annotations <- left_join(annotations, 
                         donorSamplesheet[,c("Donor", "Sex", "Predicted_Ancestry", "Age")],
                         by = "Donor") |> 
  dplyr::select(-Donor) |> 
  dplyr::rename(Ancestry = Predicted_Ancestry) |> 
  dplyr::select(Age:Condition)
annotations$Age <- cut(annotations$Age, breaks = seq(30, 90, 10),
                       labels = c("31-40", "41-50", "51-60", "61-70", "71-80", "81-90"))

# Significant genes to highlight
downGenes <- sig_degenes |>  
  filter(symbol %in% c("GREM1", "DKK1", "GDF5", "GDF10", "DLX5",
                       "GPX3", "COL21A1", "FZD8", "WWP2", "ALDH1A1",
                       "NOG", "MAP2K6", "SOX6", "NKX3-2", "ERG")) |> 
  arrange(log2FoldChange)

upGenes <- sig_degenes |>  
  filter(symbol %in% c("CXCL2", "LIF", "IL6", "IL1B", "MMP13",
                       "ADAMTS4", "NFKB1", "CAMK1G", "IRAK2", "MMP1",
                       "IL17C", "CXCR4", "WNT5A", "BMP6", "COL13A1",
                       "IL11", "CRTAC1", "COL7A1", "MMP10", "CXCL1")) |> 
  arrange(desc(log2FoldChange))

annotationObjects <- HeatmapAnnotation(
  df = annotations,
  col = list(Age = c("31-40" = ageColors[6],
                     "41-50" = ageColors[5],
                     "51-60" = ageColors[4],
                     "61-70" = ageColors[3],
                     "71-80" = ageColors[2],
                     "81-90" = ageColors[1]),
             Condition = conditionColors,
             Sex = sexColors,
             Ancestry = c("AFR" = ancestryColors[1],
                      "AMR" = ancestryColors[2],
                      "EAS" = ancestryColors[3],
                      "EUR" = ancestryColors[4],
                      "SAS" = ancestryColors[5])),
  annotation_name_gp = gpar(fontfamily = "Helvetica",
                            fontsize = 8),
  which = "column"
)

cluster_annotations <- sig_degenes |>  
  dplyr::select(gene_id, log2FC_dir) |> 
  column_to_rownames(var = "gene_id")

clusters <- HeatmapAnnotation(
  df = cluster_annotations,
  col = list(log2FC_dir = log2fcColors),
  annotation_name_gp = gpar(fontfamily = "Helvetica",
                            fontface = "bold",
                            fontsize = 0),
  simple_anno_size = unit(3, "mm"),
  which = "row"
)

# Initialize clustering of heatmap
h1 <- draw(Heatmap(mat_scaled,
                   show_row_names = FALSE,
                   top_annotation = annotationObjects,
                   right_annotation = clusters,
                   row_title = NULL,
                   cluster_columns = TRUE,
                   cluster_rows = FALSE,
                   show_column_names = FALSE,
                   col = colorRamp2(seq(-3, 3), heatmapColors),
                   show_row_dend = FALSE,
                   show_column_dend = FALSE))
# Get column order of clusters and swap CTL and FNF clusters
col_order <- column_order(h1)
new_col_order <- c(col_order[102:202], col_order[1:101])

# Plot heatmap with column order defined from clustering and CTL and FNF order above
h1 <- Heatmap(mat_scaled,
              show_row_names = FALSE,
              top_annotation = annotationObjects,
              right_annotation = clusters,
              row_title = NULL,
              cluster_columns = FALSE,
              cluster_rows = FALSE,
              show_column_names = FALSE,
              col = colorRamp2(seq(-3, 3), heatmapColors),
              show_row_dend = FALSE,
              show_column_dend = FALSE,
              column_order = new_col_order)

heatmapLegend <- Legend(at = c(-3, 3),
                        col_fun = colorRamp2(breaks = seq(-3, 3),
                                             colors = heatmapColors),
                        border = NA,
                        title_gp = gpar(fontsize = 0),
                        labels_gp = gpar(fontfamily = "Helvetica",
                                         fontsize = 8),
                        legend_width = unit(4.325, "in"),
                        grid_height = unit(0.11, "in"),
                        direction = "horizontal"
)

heatmapGrob <- grid.grabExpr(draw(h1,
                                  show_annotation_legend = FALSE,
                                  show_heatmap_legend = FALSE,
                                  background = "transparent",
                                  use_raster = TRUE))
heatmapLegendGrob <- grid.grabExpr(draw(heatmapLegend))

save(heatmapGrob, file = "plots/conditionDE_Fig1/heatmapGrob.rda")
save(heatmapLegendGrob, file = "plots/conditionDE_Fig1/heatmapLegendGrob.rda")

# GO AND KEGG BARPLOTS ----------------------------------------------------
#### GO

# Get reduced, significant GO terms for each category
upsig_go_data <- 
  read_delim("data/homer/homer_upsig_deGenes_pval01_l2fc2/biological_process.txt") |> 
  mutate(pval = exp(1)^logP) |> 
  filter(pval < 0.01)
upsig_go <- reduceGO(upsig_go_data,
                     category = "Upregulated")

## Format and write to table
upgo_table <- upsig_go |> 
  dplyr::select(-`Entrez Gene IDs`, -pval, -logP, -size, -termUniqueness, -termUniquenessWithinCluster,
                -termDispensability, -category) |> 
  relocate(`-log10pval`, .after = Enrichment) |> 
  arrange(desc(`-log10pval`))

write_csv(upgo_table, file = "tables/SupTable2A.csv")


downsig_go_data <- 
  read_delim("data/homer/homer_downsig_deGenes_pval01_l2fc2/biological_process.txt") |> 
  mutate(pval = exp(1)^logP) |> 
  filter(pval < 0.01)
downsig_go <- reduceGO(downsig_go_data,
                       category = "Downregulated")

## Format and write to table
downgo_table <- downsig_go |> 
  dplyr::select(-`Entrez Gene IDs`, -pval, -logP, -size, -termUniqueness, -termUniquenessWithinCluster,
                -termDispensability, -category) |> 
  relocate(`-log10pval`, .after = Enrichment) |> 
  arrange(desc(`-log10pval`))

write_csv(downgo_table, file = "tables/SupTable2B.csv")

# Select 5 each for plotting
upsig_go_plotting <- upsig_go |> 
  filter(Term == parentTerm) |> 
  filter(parentTerm %in% c("response to cytokine", 
                           "cell surface receptor signaling pathway", 
                           "collagen catabolic process", 
                           "regulation of cell-cell adhesion", 
                           "acute inflammatory response")) |> 
  arrange(`-log10pval`)

downsig_go_plotting <- downsig_go |> 
  filter(Term == parentTerm) |> 
  filter(parentTerm %in% c("anatomical structure development", 
                           "skeletal system development",
                           "lipid metabolic process", 
                           "positive regulation of collagen biosynthetic process",
                           "positive regulation of wound healing")) |> 
  arrange(`-log10pval`)

# Combine into one
go_plotting <- bind_rows(upsig_go_plotting, downsig_go_plotting)
go_plotting$parentTerm <- factor(go_plotting$parentTerm, levels = go_plotting$parentTerm) 
go_plotting$category <- factor(go_plotting$category, levels = c("Upregulated", "Downregulated"))

# Plot all in barplot
GO_barplots <- ggplot(go_plotting, aes(x = `-log10pval`, y = parentTerm, fill = category)) +
  geom_vline(xintercept = 10, color = "grey75", alpha = 0.4) +
  geom_vline(xintercept = 20, color = "grey75", alpha = 0.4) +
  geom_vline(xintercept = 30, color = "grey75", alpha = 0.4) +
  geom_vline(xintercept = 40, color = "grey75", alpha = 0.4) +
  geom_vline(xintercept = 50, color = "grey75", alpha = 0.4) +
  geom_bar(stat = "identity") +
  scale_x_continuous(limits = c(0, 51), expand = c(0, 0), name = "-log~10~pval",
                     breaks = seq(0, 50, 10)) +
  scale_fill_manual(values = c(log2fcColors[["+"]], log2fcColors[["-"]])) +
  facet_wrap(~category, ncol = 1, strip.position = "left", scales = "free_y") +
  geom_text(aes(x = 0, label = parentTerm), hjust = 0, family = "Helvetica",
            size = 2) +
  theme(panel.background = element_rect(fill = 'transparent', color = "transparent"),
        plot.background = element_rect(fill = 'transparent', color = "transparent"),
        text = element_text(family = "Helvetica"),
        legend.position = "None",
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_markdown(size = 6),
        axis.text.x = element_text(color = "black", size = 4),
        strip.background = element_blank(),
        axis.ticks = element_blank(),
        axis.line.x = element_line(linewidth = 0.25),
        strip.text = element_text(size = 8, color = "black"),
        panel.spacing = unit(0, "mm"), 
        plot.title = element_text(hjust = 0.5, face = "bold", size = 8)) +
  ggtitle("GO Terms")

ggsave(filename = "plots/conditionDE_Fig1/GO_barplots.pdf",
       plot = GO_barplots, width = 5, height = 8, units = "in")
save(GO_barplots, file = "plots/conditionDE_Fig1/GO_barplots.rda")

#### KEGG

# Read in from Homer
upsig_kegg_data <- read_delim("data/homer/homer_upsig_deGenes_pval01_l2fc2/kegg.txt") |> 
  mutate(pval = exp(1)^logP) |> 
  filter(pval < 0.01) |> 
  distinct(Term, .keep_all = TRUE) |> 
  mutate(`-log10pval` = -log10(pval)) |> 
  mutate(category = "Upregulated")

## Format and write to table
upsig_kegg_table <- upsig_kegg_data |> 
  dplyr::select(-logP, -pval, -`Entrez Gene IDs`, -category) |> 
  relocate(`-log10pval`, .after = Enrichment) |> 
  arrange(desc(`-log10pval`))

write_csv(upsig_kegg_table, file = "tables/SupTable2C.csv")

downsig_kegg_data <- read_delim("data/homer/homer_downsig_deGenes_pval01_l2fc2/kegg.txt") |> 
  mutate(pval = exp(1)^logP) |> 
  filter(pval < 0.01) |> 
  distinct(Term, .keep_all = TRUE) |> 
  mutate(`-log10pval` = -log10(pval)) |> 
  mutate(category = "Downregulated")

## Format and write to table

downsig_kegg_table <- downsig_kegg_data |> 
  dplyr::select(-logP, -pval, -`Entrez Gene IDs`, -category) |> 
  relocate(`-log10pval`, .after = Enrichment) |> 
  arrange(desc(`-log10pval`))

write_csv(downsig_kegg_table, file = "tables/SupTable2D.csv")

# Plot top 5 significant for each category
upsig_kegg_plotting <- upsig_kegg_data |> 
  filter(Term %in% c("TNF signaling pathway", "IL-17 signaling pathway", 
                     "Cytokine-cytokine receptor interaction", 
                     "NF-kappa B signaling pathway", "NOD-like receptor signaling pathway")) |> 
  arrange(`-log10pval`)

down_kegg_plotting <- downsig_kegg_data |> 
  filter(Term %in% c("Rap1 signaling pathway", "Drug metabolism - cytochrome P450",
                     "Calcium signaling pathway", "PPAR signaling pathway",
                     "cAMP signaling pathway")) |> 
  arrange(`-log10pval`)

kegg_plotting <- bind_rows(upsig_kegg_plotting, down_kegg_plotting)
kegg_plotting$Term <- factor(kegg_plotting$Term, levels = kegg_plotting$Term) 
kegg_plotting$category <- factor(kegg_plotting$category, levels = c("Upregulated", "Downregulated"))

KEGG_barplots <- ggplot(kegg_plotting, aes(x = `-log10pval`, y = Term, fill = category)) +
  geom_vline(xintercept = 5, color = "grey75", alpha = 0.4) +
  geom_vline(xintercept = 10, color = "grey75", alpha = 0.4) +
  geom_vline(xintercept = 15, color = "grey75", alpha = 0.4) +
  geom_vline(xintercept = 20, color = "grey75", alpha = 0.4) +
  geom_vline(xintercept = 25, color = "grey75", alpha = 0.4) +
  geom_bar(stat = "identity") +
  scale_x_continuous(expand = c(0, 0), name = "-log~10~pval", limits = c(0, 26),
                     breaks = seq(0, 25, 5)) +
  scale_fill_manual(values = c(log2fcColors[["+"]], log2fcColors[["-"]])) +
  facet_wrap(~category, ncol = 1, strip.position = "left", scales = "free_y") +
  geom_text(aes(x = 0, label = Term), hjust = 0, family = "Helvetica",
            size = 2) +
  theme(panel.background = element_rect(fill = 'transparent', color = "transparent"),
        plot.background = element_rect(fill = 'transparent', color = "transparent"),
        text = element_text(family = "Helvetica"),
        legend.position = "None",
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_markdown(size = 6),
        axis.text.x = element_text(color = "black", size = 4),
        strip.background = element_blank(),
        axis.ticks = element_blank(),
        axis.line.x = element_line(linewidth = 0.25),
        strip.text = element_blank(),
        panel.spacing = unit(0, "mm"), 
        plot.title = element_text(hjust = 0.5, face = "bold", size = 8)) +
  ggtitle("KEGG Pathways")

save(KEGG_barplots, file = "plots/conditionDE_Fig1/KEGG_barplots.rda")

# TF MOTIFS AND TF GENE EXPRESSION ----------------------------------------

de_genes_results <- read_csv("data/condition_de/de_genes_results.csv",
                             col_select = c("symbol", "log2FoldChange"))
load("data/condition_de/differential_expression_dds.rda")
normCounts <- counts(dds, normalized = TRUE)

### Upregulated motifs

# Read in and subset
upsig_knownmotifs <- read_delim("data/homer/homer_upsig_deGenes_pval01_l2fc2/knownResults.txt") |> 
  # Convert percentages to numbers
  mutate(across(c(`% of Target Sequences with Motif`, 
                  `% of Background Sequences with Motif`),
                ~ gsub("%", "", .))) |> 
  mutate(across(c(`% of Target Sequences with Motif`, 
                  `% of Background Sequences with Motif`), as.numeric)) |> 
  # Calculate log2 enrichment
  mutate(log2enrichment = log2(`% of Target Sequences with Motif`/`% of Background Sequences with Motif`)) |> 
  # Calculate -log10pval
  mutate(log10pval = -log10(exp(`Log P-value`))) |> 
  slice_max(order_by = log10pval, n = 4) |> 
  mutate(motifLogo = paste0("data/homer/homer_upsig_deGenes_pval01_l2fc2/knownResults/known", row_number(), ".logo.svg"))

# Pull out first part of motif name
upsig_knownmotifs$Name <- unlist(lapply(str_split(upsig_knownmotifs$`Motif Name`, 
                                                  '[(]'), `[[`, 1))

# Remove second NFKB
upsig_knownmotifs <- upsig_knownmotifs[-2,]

# Based on motifs, assign genes (and family members) that encode the TF
genes <- bind_rows(expand_grid(Name = upsig_knownmotifs$Name[1], 
                               Gene = c("NFKB1", "NFKB2", "RELA", "RELB", "REL")),
                   expand_grid(Name = upsig_knownmotifs$Name[2],
                               Gene = c("FOSL2", "FOS", "FOSB", "FOSL1")),
                   expand_grid(Name = upsig_knownmotifs$Name[3],
                               Gene = c("JUN", "JUNB", "JUND")))

upsig_knownmotifs <- left_join(upsig_knownmotifs, 
                               genes, 
                               by = "Name", 
                               multiple = "all")

# Get ENSEMBL IDs of genes
gene_symbols <- AnnotationDbi::select(org.Hs.eg.db, 
                                      keys = upsig_knownmotifs$Gene, 
                                      keytype = "SYMBOL", 
                                      columns = c("SYMBOL", "ENSEMBL"))

upsig_knownmotifs <- left_join(upsig_knownmotifs,
                               gene_symbols,
                               by = join_by(Gene == SYMBOL))

# For each gene, calculate each donor's FNF/CTL log2FC
upsig_donor_l2fcs <- bind_rows(lapply(unique(upsig_knownmotifs$ENSEMBL), 
                                      get_sample_l2fc, 
                                      countMatrix = normCounts))

# Join with motif data
upsig_knownmotifs_l2fc  <- left_join(upsig_knownmotifs, 
                                     upsig_donor_l2fcs, by = "ENSEMBL") |> 
  mutate(category = "up")

### Downregulated motifs

# Read in and subset for top 3
downsig_knownmotifs <- read_delim("data/homer/homer_downsig_deGenes_pval01_l2fc2/knownResults.txt") |> 
  # Convert percentages to numbers
  mutate(across(c(`% of Target Sequences with Motif`, 
                  `% of Background Sequences with Motif`),
                ~ gsub("%", "", .))) |> 
  mutate(across(c(`% of Target Sequences with Motif`, 
                  `% of Background Sequences with Motif`), as.numeric)) |> 
  # Calculate log2 enrichment
  mutate(log2enrichment = log2(`% of Target Sequences with Motif`/`% of Background Sequences with Motif`)) |> 
  # Calculate -log10pval
  mutate(log10pval = -log10(exp(`Log P-value`))) |> 
  slice_max(order_by = log10pval, n = 3) |> 
  mutate(motifLogo = paste0("data/homer/homer_downsig_deGenes_pval01_l2fc2/knownResults/known", row_number(), ".logo.svg"))

# Pull out first part of motif name
downsig_knownmotifs$Name <- unlist(lapply(str_split(downsig_knownmotifs$`Motif Name`, 
                                                    '[(]'), `[[`, 1))

# Assign gene names
downsig_knownmotifs <- downsig_knownmotifs |>  
  filter(Name == "Mef2c") |> 
  mutate(Gene = toupper(Name))

# Get ENSEMBL IDs of genes
gene_symbols <- AnnotationDbi::select(org.Hs.eg.db, 
                                      keys = downsig_knownmotifs$Gene, 
                                      keytype = "SYMBOL", 
                                      columns = c("SYMBOL", "ENSEMBL"))

downsig_knownmotifs <- left_join(downsig_knownmotifs,
                                 gene_symbols,
                                 by = join_by(Gene == SYMBOL))

# For each gene, calculate each donor's FNF/CTL log2FC
downsig_donor_l2fcs <- bind_rows(lapply(unique(downsig_knownmotifs$ENSEMBL), 
                                        get_sample_l2fc, 
                                        countMatrix = normCounts))

# Join with motif data
downsig_knownmotifs_l2fc  <- left_join(downsig_knownmotifs, 
                                       downsig_donor_l2fcs, by = "ENSEMBL") |> 
  mutate(category = "down")

# Combine data into one 
sig_knownmotifs_l2fc <- bind_rows(upsig_knownmotifs_l2fc,
                                  downsig_knownmotifs_l2fc)
sig_knownmotifs_l2fc$Name <- factor(sig_knownmotifs_l2fc$Name, 
                                    levels = c("NFkB-p65-Rel", "Fos", "JunB", "Mef2c"))
sig_knownmotifs_l2fc$category <- factor(sig_knownmotifs_l2fc$category, levels = c("up", "down"))
sig_knownmotifs_l2fc$Gene <- factor(sig_knownmotifs_l2fc$Gene, 
                                    levels = c("MEF2C", "JUND", "JUNB",
                                               "JUN", "FOS", "FOSL2", "FOSB",
                                               "FOSL1", "REL", "RELA", "NFKB1",
                                               "RELB", "NFKB2"))
tf_gene_plot <- ggplot(sig_knownmotifs_l2fc, aes(x = log2FC, y = Gene, color = category)) +
  geom_vline(xintercept = 4, color = "grey90") +
  geom_vline(xintercept = 3, color = "grey90") +
  geom_vline(xintercept = 2, color = "grey90") +
  geom_vline(xintercept = 1, color = "grey90") +
  geom_vline(xintercept = -1, color = "grey90") +
  geom_vline(xintercept = -2, color = "grey90") +
  geom_vline(xintercept = -3, color = "grey90") +
  geom_vline(xintercept = -4, color = "grey90") +
  geom_violin(color = NA, fill = "grey25", alpha = 0.2) +
  geom_jitter(size = 0.75) +
  geom_violin(color = "grey25", fill = NA, linewidth = 0.25) +
  stat_summary(fun = "median", geom = "crossbar", width = 0.5, color = "grey25", linewidth = 0.5) +
  stat_boxplot(geom = "errorbar", width = 0.25, linewidth = 0.5, color = "grey25") +
  geom_hline(data = tibble(Name = factor(c("NFkB-p65-Rel"), 
                                         levels = c("NFkB-p65-Rel", "Fos", "JunB", "Mef2c") )), aes(yintercept = Inf), linewidth = 0.25) +
  geom_hline(data = tibble(Name = factor(c("Mef2c", "JunB", "Fos"), 
                                         levels = c("NFkB-p65-Rel", "Fos", "JunB", "Mef2c") )), aes(yintercept = Inf),
             color = "grey") +
  geom_vline(data = tibble(Name = factor(c("NFkB-p65-Rel", "Fos", "JunB", "Mef2c")),
                           levels = c("NFkB-p65-Rel", "Fos", "JunB", "Mef2c")), aes(xintercept = Inf),
             linewidth = 0.25) +
  geom_vline(xintercept = 0, lty = 2, linewidth = 0.3) +
  
  scale_color_manual(values = c(log2fcColors[["+"]], log2fcColors[["-"]])) +
  scale_x_continuous(limits = c(-5, 5), expand = c(0, 0), 
                     breaks = c(-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5),
                     name = "log~2~(fold change)") +
  coord_cartesian(clip = "off") +
  facet_wrap(~Name, ncol = 1, scales = "free_y") +
  theme(panel.background = element_rect(fill = 'transparent', color = "transparent"),
        plot.background = element_rect(fill = 'transparent', color = "transparent"),
        text = element_text(family = "Helvetica"),
        legend.position = "None",
        axis.ticks = element_blank(),
        panel.spacing = unit(0, "mm"),
        axis.title.y = element_blank(),
        axis.title.x = element_markdown(size = 8),
        strip.text = element_blank(),
        axis.line = element_line(linewidth = 0.25),
        axis.text.x = element_text(color = "black", size = 6),
        axis.text.y = element_text(color = "black", size = 8))

# Convert to gtable 
tf_gene_plot_table <- ggplot_gtable(ggplot_build(tf_gene_plot))

# Resize panels to reflect number of genes in each category
panel1 <- tf_gene_plot_table$layout$t[grep('panel-1-1', tf_gene_plot_table$layout$name)]
panel2 <- tf_gene_plot_table$layout$t[grep('panel-1-2', tf_gene_plot_table$layout$name)]
panel3 <- tf_gene_plot_table$layout$t[grep('panel-1-3', tf_gene_plot_table$layout$name)]
panel4 <- tf_gene_plot_table$layout$t[grep('panel-1-4', tf_gene_plot_table$layout$name)]

tf_gene_plot_table$heights[panel1] <- 1.25 * tf_gene_plot_table$heights[panel1]
tf_gene_plot_table$heights[panel2] <- 1 * tf_gene_plot_table$heights[panel2]
tf_gene_plot_table$heights[panel3] <- 0.75 * tf_gene_plot_table$heights[panel3]
tf_gene_plot_table$heights[panel4] <- 0.25 * tf_gene_plot_table$heights[panel4]

save(tf_gene_plot_table, file = "plots/conditionDE_Fig1/tf_gene_plot_table.rda")

motifImages <- list()
# Make motif images with name and pvalues
for (motif in unique(sig_knownmotifs_l2fc$Name)){
  
  # Get image
  motifImg <- pictureGrob(readPicture(rawToChar(rsvg::rsvg_svg(sig_knownmotifs_l2fc |> 
                                                                 filter(Name == motif) |> 
                                                                 pull(motifLogo) |> 
                                                                 unique()))))
  
  motifGrab <- grid.grabExpr(expr = {
    grid.newpage()
    grid.draw(motifImg)
  })
  
  motifImages[[motif]] <- motifGrab
  
}

save(motifImages, file = "plots/conditionDE_Fig1/motifImages.rda")

# COMPARISON WITH OA EXPRESSION FROM OTHER STUDIES  ----------------------------

raak_de_genes <- read_csv("data/RAAK/RAAK_TableS3.csv", skip = 1) |> 
  dplyr::rename(symbol = GeneSYMBOL) |> 
  mutate(log2FoldChange = log2(FC),
         dir = ifelse(log2FoldChange > 0, "up", "down")) |> 
  filter(Pval < 0.05) |> 
  dplyr::select(-Pval)

fisch2018_de_genes <- read_csv("data/GSE114007/1-s2.0-S1063458418313876-mmc1.csv",
                               col_names = c("symbol", "log2FoldChange", "padj"),
                               col_select = c(1,2, 3), n_max = 12475) |> 
  filter(padj < 0.05) |> 
  mutate(dir = ifelse(log2FoldChange > 0, "up", "down")) |> 
  dplyr::select(-padj)

fu2021_de_genes <- read_csv("../DE/data/GSE168505/GSE168505_deseq_res.csv",
                            col_select = c("symbol", "log2FoldChange", "padj")) |> 
  filter(padj < 0.05) |> 
  mutate(dir = ifelse(log2FoldChange > 0, "up", "down")) |> 
  dplyr::rename(log2FoldChange_Fu = log2FoldChange,
                dir_Fu = dir) |> 
  dplyr::select(-padj)

intersection_raak_fisch_fu <- full_join(raak_de_genes, fisch2018_de_genes, by = "symbol", 
                                        suffix = c("_RAAK", "_Fisch")) |> 
  full_join(fu2021_de_genes, by = "symbol") |> 
  na.omit() |> 
  rowwise() |> 
  mutate(agree = ifelse(dir_RAAK == dir_Fisch & 
                          dir_Fisch == dir_Fu & 
                          dir_RAAK == dir_Fu, "yes", "no")) |> 
  ungroup() |> 
  filter(agree == "yes")

up_intersection_raak_fisch_fu <- intersection_raak_fisch_fu |> 
  filter(dir_RAAK == "up")

down_intersection_raak_fisch_fu <- intersection_raak_fisch_fu |> 
  filter(dir_RAAK == "down")

fnf_genes <- read_csv("data/condition_de/de_genes_results.csv")

fnf_up_oa_subset <- fnf_genes |>
  filter(symbol %in% up_intersection_raak_fisch_fu$symbol)

fnf_down_oa_subset <- fnf_genes |>
  filter(symbol %in% down_intersection_raak_fisch_fu$symbol)


fnf_all_oa_subset <- bind_rows(fnf_up_oa_subset |> mutate(group = "Up in OA"),
                               fnf_down_oa_subset |> mutate(group = "Down in OA")) |> 
  mutate(group = factor(group, levels = c("Up in OA", "Down in OA")))


up_test_oa <- wilcox.test(x = fnf_up_oa_subset$log2FoldChange,
                          y = fnf_genes |> 
                            filter(!symbol %in% fnf_up_oa_subset$symbol) |> 
                            pull(log2FoldChange),
                          alternative = "greater")

down_test_oa <- wilcox.test(x = fnf_down_oa_subset$log2FoldChange,
                            y = fnf_genes |> 
                              filter(!symbol %in% fnf_down_oa_subset$symbol) |> 
                              pull(log2FoldChange),
                            alternative = "less")

fnf_all_oa_subset <- fnf_all_oa_subset |> 
  mutate(highlight = case_when(symbol %in% c("COL3A1", "NGF", "PTGES", "ECM1", "BMP1") ~ "up",
                               symbol %in% c("WWP2", "ALDH1L1", "DDIT3", "APOD", "ABCB9") ~ "down")) |> 
  mutate(highlight = factor(highlight, levels = c("up", "down")))

## Donor-level l2fc for high confidence oa genes
load("data/condition_de/differential_expression_dds.rda")
normCounts <- counts(dds, normalized = TRUE)

goi <- c("COL3A1", "NGF", "PTGES", "ECM1", "BMP1",
         "WWP2", "ALDH1L1", "DDIT3", "APOD", "ABCB9")

highconf_geneids <- unlist(lapply(goi, symbol_to_ensembl))
highconf_gene_sample_l2fc <- lapply(highconf_geneids, get_sample_l2fc, normCounts) |> 
  bind_rows() |> 
  left_join(tibble(ENSEMBL = highconf_geneids,
                   symbol = goi),
            by = "ENSEMBL") |> 
  mutate(group = ifelse(symbol %in% c("NGF", "PTGES", "COL3A1",
                                      "ECM1", "BMP1"), "Up in OA", "Down in OA")) |> 
  mutate(group = factor(group, levels = c("Up in OA", "Down in OA")),
         symbol = factor(symbol, levels = c("ECM1", "BMP1",
                                            "COL3A1", "NGF", "PTGES",
                                            "APOD", "ABCB9", "DDIT3",
                                            "ALDH1L1", "WWP2")))
  

goi_sample_l2fc <- ggplot(highconf_gene_sample_l2fc ,
       mapping = aes(x = symbol, y = log2FC, fill = group, color = group)) +
  geom_hline(yintercept = 0, lty = 2, color = "grey25", linewidth = 0.25) +
  geom_jitter(width = 0.2, color = "grey40", size = 0.25) +
  geom_boxplot(outlier.shape = NA, linewidth = 0.5, alpha = 0.8) +
  facet_wrap(vars(group), nrow = 2, strip.position = "top", scales = "free_x") + 
  scale_fill_manual(values = c(log2fcColors[["+"]], log2fcColors[["-"]])) +
  scale_color_manual(values = c(darken(log2fcColors[["+"]], 0.3), 
                                darken(log2fcColors[["-"]], 0.3))) +
  scale_y_continuous(name = "log~2~(fold change)",
                     limits = c(-7, 7), breaks = seq(-6, 6, 2)) +
  coord_cartesian(clip = "off") +
  theme(strip.placement = "outside",
        axis.line = element_line(linewidth = 0.25),
        axis.title.x = element_blank(),
        axis.ticks = element_blank(),
        axis.ticks.length.y = unit(-0.1, "cm"),
        axis.title.y = element_markdown(size = 6, family = "Helvetica",
                                        margin = margin(r = -15)),
        text = element_text(family = "Helvetica"),
        axis.text.y = element_text(color = "black", size = 6),
        axis.text.x = element_text(color = "black", size = 8, margin = margin(b = -1)),
        strip.background = element_blank(),
        strip.text.x.bottom = element_markdown(size = 8, margin = margin(t = 1)),
        panel.background = element_rect(fill = "transparent", color = "transparent"),
        plot.background = element_rect(fill = "transparent", color = "transparent"),
        panel.grid = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, family = "Helvetica",
                                  size = 10, margin = margin(b=-3)))


save(goi_sample_l2fc, file = "plots/conditionDE_Fig1/goi_sample_l2fc.rda")

# Assemble entire figure with plotgardener --------------------------------
grDevices::cairo_pdf("plots/conditionDE_Fig1/Fig1.pdf", 
                     width = 10.75, height = 9.5)
pageCreate(width = 10.75, height = 9.5, showGuides = FALSE)

## A - DE heatmap
plotText("A", x = 0.1, y = 0.1, just = c("left", "top"), fontfamily = "Helvetica",
         fontsize = 14, fontface = "bold")
plotGG(plot = heatmapGrob, x = 0.25, y = 0.2, height = 4.5, width = 5)

# Colorbar 
plotGG(plot = heatmapLegendGrob, x = 0.325, y = 4.7, 
       width = 4.325, height = 0.11)
# Colorbar title
plotText(label = "Relative Expression", fontfamily = "Helvetica",
         fontsize = 8, x = 2.45, y = 4.8, just = "top")

# Age legend
plotRect(x = 5, y = 0.375, width = unit(3, "mm"), 
         height = unit(1, "mm"), linecolor = NA, fill = ageColors[6],
         just = "left")
plotRect(x = unit(5, "in") + unit(3, "mm"), y = 0.375, width = unit(3, "mm"), 
         height = unit(1, "mm"), linecolor = NA, fill = ageColors[5],
         just = "left")
plotRect(x = unit(5, "in") + unit(3*2, "mm"), y = 0.375, width = unit(3, "mm"), 
         height = unit(1, "mm"), linecolor = NA, fill = ageColors[4],
         just = "left")
plotRect(x = unit(5, "in") + unit(3*3, "mm"), y = 0.375, width = unit(3, "mm"), 
         height = unit(1, "mm"), linecolor = NA, fill = ageColors[3],
         just = "left")
plotRect(x = unit(5, "in") + unit(3*4, "mm"), y = 0.375, width = unit(3, "mm"), 
         height = unit(1, "mm"), linecolor = NA, fill = ageColors[2],
         just = "left")
plotRect(x = unit(5, "in") + unit(3*5, "mm"), y = 0.375, width = unit(3, "mm"), 
         height = unit(1, "mm"), linecolor = NA, fill = ageColors[1],
         just = "left")
plotText(label = "31", 
         x = 5,
         y = 0.4,
         fontfamily = "Helvetica", fontsize = 4, just = c("left", "top"))
plotText(label = "41", 
         x = unit(5, "in") + unit(3, "mm"),
         y = 0.4,
         fontfamily = "Helvetica", fontsize = 4, just = c("left", "top"))
plotText(label = "51", 
         x = unit(5, "in") + unit(3*2, "mm"),
         y = 0.4,
         fontfamily = "Helvetica", fontsize = 4, just = c("left", "top"))
plotText(label = "61", 
         x = unit(5, "in") + unit(3*3, "mm"),
         y = 0.4,
         fontfamily = "Helvetica", fontsize = 4, just = c("left", "top"))
plotText(label = "71", 
         x = unit(5, "in") + unit(3*4, "mm"),
         y = 0.4,
         fontfamily = "Helvetica", fontsize = 4, just = c("left", "top"))
plotText(label = "81", 
         x = unit(5, "in") + unit(3*5, "mm"),
         y = 0.4,
         fontfamily = "Helvetica", fontsize = 4, just = c("left", "top"))

# Ancestry Legend
plotRect(x = unit(5.05, "in") + unit(2.8*6, "mm"), 
         y = 0.575, width = unit(3, "mm"), 
         height = unit(1, "mm"), linecolor = NA, fill = ancestryColors[5],
         just = "right")
plotRect(x = unit(5.05, "in") + unit(2.8*5, "mm"), 
         y = 0.575, width = unit(3, "mm"), 
         height = unit(1, "mm"), linecolor = NA, fill = ancestryColors[4],
         just = "right")
plotRect(x = unit(5.05, "in") + unit(2.8*4, "mm"), 
         y = 0.575, width = unit(3, "mm"), 
         height = unit(1, "mm"), linecolor = NA, fill = ancestryColors[3],
         just = "right")
plotRect(x = unit(5.05, "in") + unit(2.8*3, "mm"), 
         y = 0.575, width = unit(3, "mm"), 
         height = unit(1, "mm"), linecolor = NA, fill = ancestryColors[2],
         just = "right")
plotRect(x = unit(5.05, "in") + unit(2.8*2, "mm"), 
         y = 0.575, width = unit(3, "mm"), 
         height = unit(1, "mm"), linecolor = NA, fill = ancestryColors[1],
         just = "right")
plotText(label = "SAS", x = unit(5.05, "in") + unit(2.8*5.5, "mm"),
         y = 0.65,
         fontsize = 3.5, fontfamily = "Helvetica")
plotText(label = "EUR", x = unit(5.05, "in") + unit(2.8*4.5, "mm"),
         y = 0.65,
         fontsize = 3.5, fontfamily = "Helvetica")
plotText(label = "EAS", x = unit(5.05, "in") + unit(2.8*3.5, "mm"),
         y = 0.65,
         fontsize = 3.5, fontfamily = "Helvetica")
plotText(label = "AMR", x = unit(5.05, "in") + unit(2.8*2.5, "mm"),
         y = 0.65,
         fontsize = 3.5, fontfamily = "Helvetica")
plotText(label = "AFR", x = unit(5.05, "in") + unit(2.8*1.5, "mm"),
         y = 0.65,
         fontsize = 3.5, fontfamily = "Helvetica")

# Sex legend
plotRect(x = unit(5, "in") + unit(3*6, "mm"), 
         y = 0.8, width = unit(3, "mm"), 
         height = unit(1, "mm"), linecolor = NA, fill = sexColors[["M"]],
         just = "right")
plotRect(x = unit(5, "in") + unit(3*5, "mm"), 
         y = 0.8, width = unit(3, "mm"), 
         height = unit(1, "mm"), linecolor = NA, fill = sexColors[["F"]],
         just = "right")
plotText(label = "M", x = unit(5, "in") + unit(3*5.5, "mm"),
         y = 0.875,
         fontsize = 5, fontfamily = "Helvetica")
plotText(label = "F", x = unit(5, "in") + unit(3*4.5, "mm"),
         y = 0.875,
         fontsize = 5, fontfamily = "Helvetica")

# Condition legend
plotRect(x = unit(5, "in") + unit(3*6, "mm"), 
         y = 1, width = unit(3, "mm"), 
         height = unit(1, "mm"), linecolor = NA, fill = "#4A4A4A",
         just = "right")
plotRect(x = unit(5, "in") + unit(3*5, "mm"), 
         y = 1, width = unit(3, "mm"), 
         height = unit(1, "mm"), linecolor = NA, fill = "#B8B8B8",
         just = "right")
plotText(label = "FN-f", x = unit(5, "in") + unit(3*5.5, "mm"),
         y = 1.05,
         fontsize = 4, fontfamily = "Helvetica")
plotText(label = "PBS", x = unit(5, "in") + unit(3*4.5, "mm"),
         y = 1.05,
         fontsize = 4, fontfamily = "Helvetica")

# Add labels of top downregulated and upregulated genes 
plotText(label = upGenes$symbol, x = 4.85, y = seq(1.15, 3.05, length.out = 20),
         fontcolor = "#6B5E27", fontface = "bold",
         fontsize = 6, fontfamily = "Helvetica", just = c("left", "top"))

plotText(label = downGenes$symbol, x = 4.85, y = seq(3.15, 4.5, length.out = 15),
         fontcolor = "#2F4864", fontface = "bold",
         fontsize = 6, fontfamily = "Helvetica", just = c("left", "top"))

## B - GO and KEGG
plotText("B", x = 6, y = 0.1, just = c("left", "top"), fontfamily = "Helvetica",
         fontsize = 14, fontface = "bold")
plotGG(GO_barplots, x = 5.9, y = 0.8, width = 2.5, height = 4.1)
plotGG(KEGG_barplots, x = 8.5, y = 0.8, width = 2.25, height = 4.1)

## C - TF motifs and TF gene expression
plotText("C", x = 0.1, y = 5.15, just = c("left", "top"), fontfamily = "Helvetica",
         fontsize = 14, fontface = "bold")
plotGG(tf_gene_plot_table, x = 1.575, y = 5.2, width = 3.1, height = 4.35)

plotGG(motifImages$`NFkB-p65-Rel`, x = 1.55, y = 4.8, 
       width = 1.25, height = 1.25, just = c("right", "top"))
plotText("NF-\u03BAB", x = 0.9, y = 5.575, just = "top",
         fontfamily = "Helvetica", fontsize = 8, fontface = "bold")
plotText("Upregulated", x = 0.9, y = 5.7, just = "top",
         fontsize = 8, fontfamily = "Helvetica")
plotText(paste0("pval = ", sig_knownmotifs_l2fc %>%
                  filter(Name == "NFkB-p65-Rel") %>%
                  pull(`P-value`) %>%
                  unique()),
         x = 0.9, y = 5.825, just = "top", 
         fontsize = 8, fontfamily = "Helvetica")

plotGG(motifImages$Fos, x = 1.55, y = 6.25, 
       width = 1.25, height = 1.25, just = c("right", "top"))
plotText("Fos", x = 0.9, y = 7, just = "top",
         fontfamily = "Helvetica", fontsize = 8, fontface = "bold")
plotText("Upregulated", x = 0.9, y = 7.1, just = "top",
         fontsize = 8, fontfamily = "Helvetica")
plotText(paste0("pval = ", sig_knownmotifs_l2fc %>%
                  filter(Name == "Fos") %>%
                  pull(`P-value`) %>%
                  unique()),
         x = 0.9, y = 7.215, just = "top",
         fontsize = 8, fontfamily = "Helvetica")

plotGG(motifImages$JunB, x = 1.55, 
       y = 7.45, width = 1.25, height = 1.25, just = c("right", "top"))
plotText("JunB", x = 0.9, y = 8.225, just = "top",
         fontfamily = "Helvetica", fontsize = 8, fontface = "bold")
plotText("Upregulated", x = 0.9, y = 8.318, just = "top",
         fontsize = 8, fontfamily = "Helvetica")
plotText(paste0("pval = ", sig_knownmotifs_l2fc %>%
                  filter(Name == "JunB") %>%
                  pull(`P-value`) %>%
                  unique()),
         x = 0.9, y = 8.435, just = "top",
         fontsize = 8, fontfamily = "Helvetica")

plotGG(motifImages$Mef2c, x = 1.55, y = 8.3, 
       width = 1.25, height = 1.25, just = c("right", "top"))
plotText("Mef2c", x = 0.9, y = 9.05, just = "top",
         fontfamily = "Helvetica", fontsize = 8, fontface = "bold")
plotText("Downregulated", x = 0.9, y = 9.148, just = "top",
         fontsize = 8, fontfamily = "Helvetica")
plotText(paste0("pval = ", sig_knownmotifs_l2fc %>%
                  filter(Name == "Mef2c") %>%
                  pull(`P-value`) %>%
                  unique()),
         x = 0.9, y = 9.265, just = "top",
         fontsize = 8, fontfamily = "Helvetica")

## D - Comparison to OA
plotText("D", x = 5, y = 5.15, just = c("left", "top"), fontfamily = "Helvetica",
         fontsize = 14, fontface = "bold")
oa_boxplot_plot(data = fnf_all_oa_subset, up_wilcox_test = up_test_oa,
                down_wilcox_test = down_test_oa, 
                x = unit(5.5, "native"), 
                y = unit(5.6, "native"), width = 3.75, height = 3)
plotGG(goi_sample_l2fc, x = 7.85, y = 5.15, width = 2.75, height = 4.4)

dev.off()

