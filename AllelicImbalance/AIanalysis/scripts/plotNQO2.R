library(plotgardener)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)


pdf(file = "plots/NQO2_AI_snps.pdf", width = 6, height = 2.5)
pageCreate(width = 6, height = 2.5, showGuides = FALSE)

genePlot <- plotGenes(params = pgParams(gene = "NQO2", assembly = "hg38"), 
                      x = 0, y = 0.5, width = 6, height = 1.5, 
                      geneHighlights = data.frame("gene" = c("NQO2"), "color" = c("black")))


rsCTL <- pgParams(chrom = "chr6", chromstart = 3016904, chromend = 3016904)
rsFNF <- pgParams(chrom = "chr6", chromstart = 3018703, chromend = 3018703)

annoHighlight(plot = genePlot, params = rsCTL, y = .85, height = .425,
              fill = "#48A9DF", linecolor = "#48A9DF", lwd = 1.5, alpha = 1)
annoHighlight(plot = genePlot, params = rsFNF, y = .85, height = .425,
              fill = "#EF734A", linecolor = "#EF734A", lwd = 1.5, alpha = 1)
dev.off()
