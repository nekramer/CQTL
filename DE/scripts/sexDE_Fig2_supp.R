library(tidyverse)
source("../plotting_utils.R")

# Plotting function -------------------------------------------------------

plotSexResponseGene <- function(geneRow, dds){
  
  # Gene ID
  geneid <- geneRow[6]
  # Gene symbol
  genesymbol <- geneRow[13]
  # MvsF Contrast stats
  padj <- format(as.numeric(geneRow[12]), digits = 3, scientific = TRUE)
  lfc <- round(as.numeric(geneRow[8]), digits = 3)
  
  # F CTLvsFNF stats
  F_stats <- as.data.frame(results(dds, name = "SexF.ConditionFNF")) %>% 
    rownames_to_column(var = "gene_id") %>%
    filter(gene_id == geneid)
  F_padj <- format(F_stats %>% pull(padj), digits = 3, scientific = TRUE)
  F_lfc <- round(F_stats %>% pull(log2FoldChange), digits = 3)
  
  F_df <- data.frame("Sex" = "F",
                     "count" = 100)
  
  # M CTLvsFNF stats
  M_stats <- as.data.frame(results(dds, name = "SexM.ConditionFNF")) %>% 
    rownames_to_column(var = "gene_id") %>%
    filter(gene_id == geneid)
  M_padj <- format(M_stats %>% pull(padj), digits = 3, scientific = TRUE)
  M_lfc <- round(M_stats %>% pull(log2FoldChange), digits = 3)
  M_df <- data.frame("Sex" = "M",
                     "count" = 100)
  
  
  # Grab count data for gene from plotCounts
  geneData <- plotCounts(dds = dds, gene = geneid, intgroup = c("Sex", "Condition"),
                         returnData = TRUE) %>%
    mutate(group = paste0(Sex, "_", Condition))
  
  
  # Plot with ggplot
  plot <- ggplot(geneData, aes(x = Condition, y = count, color = Sex)) +
    geom_jitter(width = 0.1, size = 2) +
    facet_wrap(~Sex, strip.position = "bottom") +
    scale_y_continuous(trans = "log2", breaks = c(100, 200, 500, 1000, 2000), 
                       name = "Normalized counts") +
    scale_x_discrete(labels = c("CTL", "FNF", "CTL", "FNF")) +
    scale_color_manual(values = c(sexColors[["F"]], sexColors[["M"]])) +
    geom_richtext(data = F_df, aes(x = 1.5, family = "Helvetica"),
                  size = 3,
                  label = paste0("padj = ", F_padj, "<br>", "log~2~FC = ", F_lfc),
                  fill = NA, label.color = NA) +
    geom_richtext(data = M_df, aes(x = 1.5, family = "Helvetica"),
                  size = 3,
                  label = paste0("padj = ", M_padj, "<br>", "log~2~FC = ", M_lfc),
                  fill = NA, label.color = NA) +
    theme_custom_scatterplot() +
    theme(legend.position = "none",
          axis.title.x = element_blank(),
          strip.placement = "outside",
          strip.background = element_blank(),
          strip.text = element_text(size = 10),
          strip.switch.pad.wrap = unit(0, "in"),
          plot.subtitle = element_markdown(hjust = 0.05, size = 10, color = "grey25"),
          plot.title = element_text(hjust = 0.05)) +
    ggtitle(label = genesymbol, subtitle = paste0("padj = ", padj, "<br>",
                                                  "log~2~FC = ", lfc))
  
  # Save
  ggsave(plot, filename = paste0("plots/sexDE_Fig2_supp/sexspecifictreatment_",
                                 genesymbol, ".pdf"),
         width = 6, height = 6, units = "in" )
  return(plot)
  
} 

# Load and plot genes -----------------------------------------------------
load("data/sex_de/dds_sex_treatment_response.rda")
sexDE_treatmenteffect_pval01 <- 
  read_csv("data/sex_de/sexDE_treatmenteffect_pval01.csv")

apply(sexDE_treatmenteffect_pval01, 
      1, plotSexResponseGene, dds = dds_sex_treatment_response)
