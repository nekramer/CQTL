library(tidyverse)
library(DESeq2)
library(mariner)

# Load extracted loops counts ---------------------------------------------
loopCounts <- readRDS("data/hic/all_loopCounts_10kb.rds")

# Concatenate loopCounts and DESeqResults ---------------------------------
load("data/hic/hic_dds.Rdata")

hic_res <- lfcShrink(dds, coef = "Treatment_PBS_vs_FNF", 
                 type = "apeglm")
save(hic_res, file = "data/hic/hic_res.Rdata")

mcols(interactions(loopCounts)) <- assay(loopCounts)
mcols(interactions(loopCounts)) <- cbind(mcols(interactions(loopCounts)), hic_res)

diff_loop_gi <- loopCounts |> 
  interactions() |> 
  as.data.frame() |> 
  as_ginteractions() |> 
  subset(!is.na(padj))


CQTL_10kb_static_loops <- diff_loop_gi |> 
  subset(padj >= 0.1)
save(CQTL_10kb_static_loops, file = "data/hic/CQTL_10kb_static_loops.rda")

CQTL_10kb_sig_gained_loops <- diff_loop_gi |> 
  subset(padj < 0.1 & log2FoldChange < 0)
save(CQTL_10kb_sig_gained_loops, file = "data/hic/CQTL_10kb_sig_gained_loops.rda")
