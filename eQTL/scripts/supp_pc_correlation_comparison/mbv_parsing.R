#!/usr/bin/R
library(tidyverse)


# Arguments ---------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
mbv <- read_delim(args[1])
sample <- args[2]

# Parse sample name for donor
rnaDonor <- unlist(str_split(sample, "_"))[2]
condition <- unlist(str_split(sample, "_"))[4]

# Calculate concordance at hets and concordance at homs -------------------

mbv <- mbv %>% 
  mutate(hetConcordance = n_het_consistent/n_het_covered,
         homConcordance = n_hom_consistent/n_hom_covered)


# Get matching donor based on maximum het and hom concordance ------------

matchingDonor <- mbv[which(mbv$hetConcordance == max(mbv$hetConcordance) & 
                             mbv$homConcordance == max(mbv$homConcordance)),] %>%
  pull(SampleID)

# Check against donor of bam file
sampleConcordance <- tibble("Donor" = rnaDonor,
                            "Condition" = condition,
                            "mbvMatch" = matchingDonor == rnaDonor) %>%
  write_csv(file = args[3])


