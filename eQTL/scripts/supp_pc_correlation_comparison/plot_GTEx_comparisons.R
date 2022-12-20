library(tidyverse)
library(RColorBrewer)

args <- commandArgs(trailingOnly = TRUE)


# plot eGene pi1 ----------------------------------------------------------
# Read in CTL, add Condition column, arrange by pi1, then convert tissue order to factor
CTL_eGene_pi1 <- read_csv(args[1]) %>% mutate(Condition = "CTL") %>%
  arrange(pi1) %>%
  # Clean tissue labels
  mutate(across(tissue, gsub, pattern = "_", replacement = " ")) %>%
  mutate(across(tissue, factor, levels = .$tissue))
# Read in FNF, add Condition column, convert tissue order to factor with CTL levels
FNF_eGene_pi1 <- read_csv(args[2]) %>% mutate(Condition = "FNF") %>%
  # Clean tissue labels
  mutate(across(tissue, gsub, pattern = "_", replacement = " ")) %>%
  mutate(across(tissue, factor, levels = CTL_eGene_pi1$tissue))
  
# Combine CTL and FNF for plotting
eGene_pi1 <- bind_rows(CTL_eGene_pi1, FNF_eGene_pi1) %>%
  # Convert Condition column to factor
  mutate(across(Condition, as.factor))
  

# Plot pi1 estimates
ggplot(eGene_pi1, mapping = aes(x = pi1, y = tissue)) +
  geom_point(aes(color = Condition)) +
  theme_minimal() +
  xlim(0, 1) +
  xlab("pi1") +
  #xlab("\u03C0\u2081") +
  scale_color_manual(values = c(brewer.pal(n = 6, "YlGnBu")[3], brewer.pal(n = 6, "YlGnBu")[5]),
                     labels = c("Control", "FN-f")) +
  theme(axis.title.y=element_blank(),
        legend.title = element_blank())

ggsave(filename = "output/GTEx/GTEx_eGene_pi1.pdf")


# downsampled eGene percent overlaps --------------------------------------

CTL_eGene_percOverlap <- read_csv(args[3]) %>% mutate(Condition = "Control") %>%
  arrange(percentOverlap) %>%
  # Clean tissue labels
  mutate(across(tissue, gsub, pattern = "_", replacement = " ")) %>%
  mutate(across(tissue, factor, levels = .$tissue))
  
  
FNF_eGene_percOverlap <- read_csv(args[4]) %>% mutate(Condition = "FN-f") %>%
  # Clean tissue labels
  mutate(across(tissue, gsub, pattern = "_", replacement = " ")) %>%
  mutate(across(tissue, factor, levels = CTL_eGene_percOverlap$tissue))

# Combine CTL and FNF for plotting
eGene_percOverlap <- bind_rows(CTL_eGene_percOverlap, FNF_eGene_percOverlap) %>%
  # Convert Condition column to factor
  mutate(across(Condition, as.factor)) %>%
  mutate(percent = percentOverlap*100)

# Plot downsampled percent Overlap
ggplot(eGene_percOverlap, mapping = aes(x = percent, y = tissue)) +
  geom_bar(stat = "identity") +
  facet_wrap(~Condition,
             scales = "free_x") +
  theme_minimal() +
  xlab("eGene Percent Overlap") +
  theme(axis.title.y=element_blank(),
        legend.title = element_blank(),
        strip.text.y = element_blank())
ggsave(filename = "output/GTEx/GTEx_eGene_percentOverlap.pdf")