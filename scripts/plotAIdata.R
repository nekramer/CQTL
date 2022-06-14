library(ggplot2)
library(DESeq2)

plotAIdata <- function(snp, gene, dds){
    
    # Grab data 
    snpData <- plotCounts(dds = dds, gene = snp,
                          intgroup = c("Treatment", "Donor", "Allele"),
                          returnData = TRUE)
    
    # Plot
    plot <- ggplot(data = snpData, aes(Allele, count, color = Treatment)) + 
        scale_color_manual(values = c("#48A9DF", "#EF734A")) +
        theme_minimal() +
        geom_hline(aes(yintercept = 0), color = "grey") +
        geom_line(aes(group = Donor), color = "grey", lwd = 0.1) +
        geom_point() + 
        facet_wrap(~Treatment, strip.position = "bottom", 
                   labeller = as_labeller(c("CTL" = "Control",
                                            "FNF" = "FN-f"))) +
        theme(panel.grid = element_blank(),
              axis.title.x = element_blank(),
              axis.text.x.bottom = element_text(vjust = 20),
              strip.text.x = element_text(vjust = -2),
              legend.position = "none",
              axis.ticks.y = element_line(color = "grey")) +
        labs(title = gene, subtitle = snp)
    return(plot)
}