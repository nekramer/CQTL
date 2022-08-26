library(tidyverse)

variantMismatches <- read_csv('reports/variantMismatches.csv')

mismatchedDonor_cumsum <- ggplot(data = variantMismatches, mapping = aes(x = `# Geno/RNA mismatched donors`)) + 
  stat_bin(aes(y=cumsum(..count..)), bins = 50, geom = "step") + 
  theme_minimal()

ggsave('plots/mismatchedDonor_cumsum.png', mismatchedDonor_cumsum)

# Number of variants removed based on having any mismatching donors
mismatching_no <- variantMismatches %>% 
  filter(`# Geno/RNA mismatched donors` > 0) %>%
  nrow()

meanMismatchedDonor <- mean(variantMismatches$`# Geno/RNA mismatched donors`)

mismatchedDonor_hist <- ggplot(data = variantMismatches, mapping = aes(x = `# Geno/RNA mismatched donors`)) +
  geom_histogram() +
  theme_minimal()
ggsave('plots/mismatchedDonor_hist.png', mismatchedDonor_hist)


genotypingErrors <- read.csv('reports/GenoRNA_errorStats.csv')
allMatchingSubset <- na.omit(genotypingErrors)
rangeGenoHet <- range(allMatchingSubset$X..genoHet_anyZeroCounts)


genoHom_10Alt <- allMatchingSubset[which(allMatchingSubset$X..genoHom_10AltAlleleCounts > 0),] %>% nrow()

remove_v1 <- read_csv('reports/remove_v1.csv')
total_removed_v1 <- remove_v1 %>% nrow()


total_removed_afterHetThreshold_v1 <- read_csv('reports/removeVariants_v1_afterHetThreshold.csv')
hetThreshold_removed <- total_removed_afterHetThreshold_v1$variant[!total_removed_afterHetThreshold_v1$variant %in% remove_v1$variant]


genotypingErrors_v2 <- read_csv('reports/GenoRNA_errorStats_v2.csv')
genoHom_10Alt_v2 <- genotypingErrors_v2[which(genotypingErrors_v2$`# genoHom_10AltAlleleCounts` > 0),] %>% nrow()

remove_v2 <- read_csv('reports/remove_v2.csv')

hetThreshold_removed_v2 <- read.csv('reports/hetThresholdRemove_v2.csv')
total_removed_afterHetThreshold_v2 <- nrow(hetThreshold_removed_v2) + nrow(remove_v2)


hetNumbers_v2 <- read_csv('reports/hetNumbers_v2.csv')

combined_stats <- variantMismatches %>% 
  left_join(genotypingErrors_v2) %>% 
  left_join(hetNumbers_v2) %>%
  left_join(rnaHets) %>%
  left_join(genoHets)

