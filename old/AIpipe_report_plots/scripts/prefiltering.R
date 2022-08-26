library(tidyverse)
snpStats <- read.table("reports/snpStats.txt", header = TRUE)

numDonors <- snpStats$group %>%
  str_split("_") %>% 
  lapply(`[[`, 2) %>% 
  unlist() %>%
  unique() %>% 
  length()

snphist <- ggplot(data = snpStats, mapping = aes(x = numSNPs)) +
  geom_histogram(bins = 40) +
  theme_minimal() +
  scale_x_continuous(breaks = scales::extended_breaks(n = 10),
                     limits = c(600000, 1350000)) +
  annotate(geom = "text", x= 1150000, y = 8, label = paste(numDonors, "donors"),
           size = 8) + 
  ylab(label = "counts")

ggsave(filename = "plots/snphist.png", plot = snphist)

maxSNPs <- max(snpStats$numSNPs)
minSNPs <- min(snpStats$numSNPs)

sample_maxSNPs <- snpStats[which(snpStats$numSNPs == maxSNPs), "group"]
sample_minSNPs <- snpStats[which(snpStats$numSNPs == minSNPs), "group"]


unfiltered_commonSnps <- read_csv("reports/unfiltered_commonSnps.txt", 
                                  col_names = FALSE) %>%
  str_remove("Number of unfiltered common SNPs: ") %>%
  as.numeric()

commonSNPs_maxDif <- maxSNPs - unfiltered_commonSnps
commonSNPs_minDif <- minSNPs - unfiltered_commonSnps


######### Assessing which variants are common
alleleCount_stats <- read.table("reports/alleleCount_stats.txt", header = TRUE)
total_unique <- nrow(alleleCount_stats)

perc_removed <- (total_unique - unfiltered_commonSnps)/total_unique * 100
ggplot(data = alleleCount_stats, mapping = aes(x = num_groups)) + 
  geom_histogram() +
  theme_minimal()


