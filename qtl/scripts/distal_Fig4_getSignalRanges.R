# This is the script used to define the signal ranges of each of the independent
# signals of significant eGenes. The signal range is defined as the variants
# in a signal with LD > 0.6 with the lead variant.
# These ranges are appended to the top signal data with the columns
# `min_signal_pos` and `max_signal_pos`.

library(tidyverse)
library(data.table)
# Define a function to do this for any condition prefix
define_signalRanges <- function(name_prefix){
  
  # Read in conditionally independent signals with their leads and corresponding LD buddies
  signals_topVariants_LD_ranges <- fread(paste0("data/eqtl/", name_prefix,
                                                "_cond1Mb_topSignals_rsID_LD.csv"),
                                  data.table = FALSE) |> 
    # Filter for LD R2 values > 0.6
    filter(R2 > 0.6) |> 
    # Extract ld buddy position from ld_variantID
    separate_wider_delim(cols = "ld_variantID", delim = ":", 
                         names = c(NA, "ld_pos", NA, NA), cols_remove = FALSE) |> 
    mutate(ld_pos = as.numeric(ld_pos)) |> 
    # Group by top variant per signal
    group_by(gene_id, variantID) |> 
    # Get minimum and maximum ld_pos
    summarize(min_signal_pos = min(ld_pos),
              max_signal_pos = max(ld_pos), 
              .groups = "keep")
  
  # Join min_signal_pos and max_signal_pos back to top variant signal data
  signals_topVariants_signalRanges <- read_csv(paste0("data/eqtl/", name_prefix,
                                         "_cond1Mb_topSignals_rsID.csv")) |> 
    left_join(signals_topVariants_LD_ranges, by = c("gene_id", "variantID"))
  
  # Write to file
  write_csv(signals_topVariants_signalRanges, file = paste0("data/eqtl/",
                                                            name_prefix, 
                                                            "_cond1Mb_topSignals_rsID_signalRanges.csv"))
}



invisible(lapply(c("CTL_PEER_k20_genoPC", "FNF_PEER_k22_genoPC"), 
                 define_signalRanges))
  

    

