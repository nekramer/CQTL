library(tidyverse)
library(data.table)

totalAlleleCountsReport <- fread("reports/totalAlleleCountsReport_10_2.csv")

# Confirming number of variants based on groupings:
totalAlleleCountsReport %>% group_by(donor, condition) %>% count() %>% ungroup()


## For each variant, do the different donor/condition combos agree in terms of 
## both thresholds?

## For the ones that agree for both thresholds, how many satisfy the min 
## total Allele count threshold? How many satisfy the min allele threshold?
## How many satisfy both?

byVariant <- totalAlleleCountsReport %>% 
  group_by(variant) %>% 
  group_split()
  
vars <- totalAlleleCountsReport %>% 
  group_by(variant) %>% 
  group_keys()
  
names(byVariant) <- vars$variant


# Define a function to count how many agree for the first and second thresholds
donorAgreement <- function(varTibble){
  summarizedDonor <- varTibble %>%
    group_by(donor) %>%
    mutate(same1 = +(n_distinct(`>= minTotalAlleleCounts`)),
           same2 = +(n_distinct(`>= minAlleleThreshold`))) %>%
    ungroup()
  
  minTotalAlleleCounts_agree <- summarizedDonor[which(summarizedDonor$same1 == 1),] %>% nrow() %>% "/"(2)
  minTotalAlleleThreshold_agree <- summarizedDonor[which(summarizedDonor$same2 == 1),] %>% nrow() %>% "/"(2)
  number_bothAgree <- summarizedDonor[which(summarizedDonor$same1 == 1 & summarizedDonor$same2 == 1),] %>% nrow() %>% "/"(2)
  return(c(minTotalAlleleCounts_agree, minTotalAlleleThreshold_agree, number_bothAgree))
}

## Ran this once since it takes a while and saved it
thresholdAgreements <- lapply(byVariant, donorAgreement)
save(thresholdAgreements, file = "reports/thresholdAgreements.rda")

load(file = "reports/thresholdAgreements.rda")

# Convert to data frame
thresholdAgreements <- as.data.frame(do.call(rbind, thresholdAgreements))
colnames(thresholdAgreements) <- c("minTotalAlleleCounts", 
                                               "minAlleleThreshold",
                                   "bothAgree")

minTotal_plot <- ggplot(data = thresholdAgreements, mapping = aes(x = minTotalAlleleCounts)) +
  geom_histogram(bins = 20) +
  theme_minimal() +
  xlab(label = "Number of donors that agree for minTotalAlleleCount logic across conditions")

ggsave(filename = "plots/minTotal_plot.png", plot = minTotal_plot)

minAllele_plot <- ggplot(data = thresholdAgreements, mapping = aes(x = minAlleleThreshold)) +
  geom_histogram() +
  theme_minimal() +
  xlab(label = "Number of donors that agree for minAlleleThreshold logic across conditions")

ggsave(filename = "plots/minAllele_plot.png", plot = minAllele_plot)

# Number of variants with all donors having agreeing minTotalAlleleCounts
var_agreeMinTotal <- length(which(thresholdAgreements$minTotalAlleleCounts == 56))
var_agreeMinTotal_perc <- var_agreeMinTotal/nrow(thresholdAgreements)*100
# Number of variants with all donors having agreeing minAlleleThresholdCounts
var_agreeAlleleTotal <- length(which(thresholdAgreements$minAlleleThreshold == 56))
var_agreeAlleleTotal_perc <- var_agreeAlleleTotal/nrow(thresholdAgreements)*100

bothThresholds_plot <- ggplot(data = thresholdAgreements, mapping = aes(x = bothAgree)) +
  geom_histogram() +
  theme_minimal() +
  scale_y_continuous(limits = c(0, 60000)) +
  xlab(label = "Number of donors that agree for both sets of logic across conditions")
ggsave(filename = "plots/bothThresholds_plot.png", plot = bothThresholds_plot)

# Define a function to get stats from agreeing donor conditions (i.e. from this pool, how many can we call as het?)
donorAgreement_subset <- function(varTibble){
  summarizedDonor <- varTibble %>%
    group_by(donor) %>%
    mutate(same1 = +(n_distinct(`>= minTotalAlleleCounts`)),
           same2 = +(n_distinct(`>= minAlleleThreshold`))) %>%
    ungroup()
  bothAgree <- summarizedDonor[which(summarizedDonor$same1 == 1 & summarizedDonor$same2 == 1),]
  summarizedDonor_subset <- distinct(bothAgree, donor, .keep_all = TRUE)
  num_agree <- summarizedDonor_subset %>% nrow()
  
  num_minTotal_het <- summarizedDonor_subset[which(summarizedDonor_subset$`>= minTotalAlleleCounts`),] %>% nrow()
  num_minAllele_het <- summarizedDonor_subset[which(summarizedDonor_subset$`>= minAlleleThreshold`),] %>% nrow()
  # nHets <- summarizedDonor_subset[which(summarizedDonor_subset$`>= minTotalAlleleCounts` & summarizedDonor_subset$`>= minAlleleThreshold`),] %>%
  #   nrow()
  return(c(num_agree, num_minTotal_het, num_minAllele_het))
  
}

## Ran this once since it takes a while and saved it
thresholdSatisfiers <- lapply(byVariant, donorAgreement)
save(thresholdSatisfiers, file = "reports/thresholdSatisfiers.rda")

load("reports/thresholdSatisfiers.rda")

rnaHets <- read.table('reports/total_RNAhets.csv', sep = ",", header = TRUE)
## Put rnaHets in same order as thresholdSatisfiers
rnaHets <- rnaHets[match(names(thresholdSatisfiers), rnaHets$variant),]

thresholdSatisfiers <- as.data.frame(do.call(rbind, thresholdSatisfiers))
colnames(thresholdSatisfiers) <- c("numAgreeingDonors", 
                                   "minTotalTrue",
                                   "minAlleleTrue")


rnaHets$percentage_valid_hets <- rnaHets$rnahet/thresholdSatisfiers$numAgreeingDonors*100


heterozygote_percentage <- ggplot(data = rnaHets, mapping=aes(x = percentage_valid_hets)) +
  geom_histogram() +
  theme_minimal() +
  scale_y_continuous(limits = c(0, 60000)) +
  xlab(label = "Percentage of heterozygotes (out of agreeing thresholds donor/conds)")
ggsave(filename = "plots/heterozygote_percentage.png", plot = heterozygote_percentage)
