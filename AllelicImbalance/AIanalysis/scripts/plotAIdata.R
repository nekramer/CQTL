library(ggplot2)
library(DESeq2)

plotAIdata <- function(snp, gene, dds, rsid){
    
    # Grab data and remove homozygous donors from plotting
    snpData <- plotCounts(dds = dds, gene = snp,
                          intgroup = c("condition", "donor", "allele"),
                          returnData = TRUE) %>%
      group_by(donor, condition) %>% 
      filter(sum(count) >= 10) %>% 
      filter(all(count >= 2)) %>%
      ungroup()
      
    # Plot
    plot <- ggplot(data = snpData, aes(allele, count, color = condition)) + 
        scale_color_manual(values = c("#48A9DF", "#EF734A")) +
        theme_minimal() +
        geom_hline(aes(yintercept = 0), color = "grey") +
        geom_line(aes(group = donor), color = "grey", lwd = 0.1) +
        geom_point() + 
        facet_wrap(~condition, strip.position = "bottom", 
                   labeller = as_labeller(c("CTL" = "Control",
                                            "FNF" = "FN-f"))) +
        ylab("normalized RNA-seq counts") +
        theme(panel.grid = element_blank(),
              axis.title.x = element_blank(),
              axis.text.x.bottom = element_text(vjust = 15),
              strip.text.x = element_text(vjust = -1, size = 12),
              legend.position = "none",
              axis.ticks.y = element_line(color = "grey"),
              plot.title = element_text(face = "bold")) +
        labs(title = gene, subtitle = rsid)
    return(plot)
}

plotDonorAltFraction <- function(snp, gene, dds, rsid){
  snpData <- plotCounts(dds = dds, gene = snp,
                        intgroup = c("condition", "donor", "allele"),
                        returnData = TRUE) %>%
    group_by(donor, condition) %>% mutate(total = sum(count)) %>% 
    filter(total >= 10) %>%
    filter(all(count >= 2)) %>%
    ungroup() %>%
    filter(allele == "alt") %>%
    mutate(alt_fraction = count/total)
  if (nrow(snpData) == 0){
    return(NULL)
  }
  
  plot <- ggplot(data = snpData, aes(x = donor, y = alt_fraction, color = condition)) + 
    scale_color_manual(values = c("#48A9DF", "#EF734A"), labels = c("Control", "FN-f")) +
    geom_point() + 
    theme_minimal() +
    ylim(0, 1) + ylab("alt allele reads/total") +
    geom_hline(aes(yintercept = 0.5), color = "grey", linetype = 2) +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 45),
          legend.title = element_blank()) +
    labs(title = gene, subtitle = rsid)
  return(plot)
  
}
