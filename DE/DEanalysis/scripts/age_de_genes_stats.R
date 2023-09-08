library(DESeq2)
library(tidyverse)
library(rrvgo)
source("scripts/plotting_utils.R")


splineFitAge <- function(gene, coef_mat, design_mat){

  gene_spline_fit <- data.frame("Sample" = rownames(design_mat),
                                "logmu" = design_mat %*% coef_mat[gene,],
                                "gene_id" = gene) %>%
    separate_wider_delim(cols = "Sample", delim = "_", 
                         names = c(NA, NA, NA, NA, NA, NA, "Age"),
                         cols_remove = FALSE) %>%
    mutate(Age = as.numeric(Age))
    
    return(gene_spline_fit)
  
}



# Union of untreated and treated genes ------------------------------------

ctl_cluster_pval01 <- read_csv("data/ctl_age_pval01clusters.csv") %>%
  distinct(gene_id, .keep_all = TRUE) %>%
  dplyr::select(gene_id, cluster)


fnf_cluster_pval01 <- read_csv("data/fnf_age_pval01clusters.csv") %>%
  distinct(gene_id, .keep_all = TRUE) %>%
  dplyr::select(gene_id, cluster)

# Join and order by cluster
union_sig_genes <- full_join(ctl_cluster_pval01, fnf_cluster_pval01) %>%
  arrange(desc(cluster))


# Get spline fit data from control ----------------------------------------

load("data/dds_age_ctl_lrt.rda")

coef_mat <- coef(dds_age_ctl_lrt)
design_mat <- model.matrix(design(dds_age_ctl_lrt), colData(dds_age_ctl_lrt))

ctl_sig_gene_fits <- bind_rows(lapply(union_sig_genes$gene_id, splineFitAge, coef_mat, design_mat)) %>%
  left_join(union_sig_genes, by = "gene_id")

# Full cluster line plots ------------------------------------------------------

for (clust in 1:2){
  
  if (clust == 1){
    median_color <- "#009E8E"
    limits <- c(0, 16)
    breaks <- seq(4, 16, 4)

  } else if (clust == 2){
    median_color <- "#AC7F21"
    limits <- c(0,14)
    breaks <- seq(2, 14, 4) 
  } 
  
  ctl_cluster <- ctl_sig_gene_fits %>%
    filter(cluster == clust)
  
  ggplot(ctl_cluster, aes(x = Age, y = logmu, color = gene_id)) +
    geom_line(lwd = 0.5) +
    scale_color_manual(values = rep("grey75", nrow(ctl_cluster))) +
    stat_summary(geom = "line", fun = "median", color = median_color, lwd = 1) +
    geom_vline(xintercept = 62, lty = 2) +
    scale_x_continuous(breaks = c(40, 50, 60, 70, 80)) +
    scale_y_continuous(name = "Log Expression Fit",
                       limits = limits,
                       breaks = breaks) +
    coord_cartesian(clip = "off") +
    theme(panel.background = element_blank(),
          text = element_text(family = "Helvetica"),
          panel.border = element_rect(fill = NA),
          legend.position = "none",
          axis.text = element_text(color = "black", size = 12),
          axis.title = element_text(size = 12),
          axis.line = element_line(color = "black"),
          axis.ticks.x = element_blank(),
          axis.ticks.length.y = unit(-0.1, "cm"))
  
  ggsave(filename = paste0("plots/splineFit_unionsigGenes_cluster", clust, "_splines.pdf"),
         width = 6, height = 6, units = "in")
  
}



