library(tidyverse)

genoHets <- read.table('reports/genohetStats.csv', header = TRUE)
rnaHets <- read.table('reports/total_RNAhets.csv', sep = ",", header = TRUE)

genoHet_hist <- ggplot(data = genoHets, mapping = aes(x = totalGenohets)) +
  geom_histogram(bins = 40) +
  theme_minimal() + 
  scale_y_continuous(limits = c(0, 35000)) +
  xlab(label = "total genotyping hets")

ggsave(filename = "plots/genoHet_hist.png", plot = genoHet_hist)


rnaHet_hist <- ggplot(data = rnaHets, mapping = aes(x = rnahet)) +
  geom_histogram(bins = 40) +
  theme_minimal() +
  scale_y_continuous(limits = c(0, 35000)) +
  xlab(label = "total RNA hets")

ggsave(filename = "plots/rnaHet_hist.png", plot = rnaHet_hist)


combinedHets <- genoHets %>% left_join(rnaHets) 

numVariants_mismatching <- combinedHets[which(combinedHets$totalGenohets != combinedHets$rnahet),] %>% nrow()
perc <- (numVariants_mismatching/nrow(combinedHets))*100

number_var_with0_genoHets <- combinedHets[which(combinedHets$totalGenohets == 0),] %>% nrow()
number_var_with0_rnaHets <- combinedHets[which(combinedHets$rnahet == 0),] %>% nrow()
