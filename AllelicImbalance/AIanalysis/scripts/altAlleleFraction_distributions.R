library(DESeq2)
library(RColorBrewer)
library(tidyr)
library(ggplot2)
library(readr)
source("scripts/utils.R")


# FUNCTIONS ---------------------------------------------------------------

get_altFractions <- function(data, dds){
  snp <- data["snp"]
  rsid <- data["rsid"]
  snpData_alt <- plotCounts(dds = dds, gene = snp,
                            intgroup = c("condition", "donor", "allele"),
                            returnData = TRUE) %>%
    group_by(donor, condition) %>% mutate(total = sum(count)) %>% 
    filter(total >= 10) %>%
    filter(all(count >= 2)) %>%
    ungroup() %>%
    filter(allele == "alt") %>%
    mutate(alt_fraction = count/total) %>%
    mutate(snp = snp) %>%
    mutate(rsid = rsid)
  # Try a mean
  mean_altAllele <- data.frame("snp" = snp,
                               "rsid" = rsid,
                               "mean_alt_Fraction" = mean(snpData_alt$alt_fraction))
  
  
  #return(snpData_alt)
  return(mean_altAllele)
}

get_condition_altFractions <- function(data, dds, cond){
  snp <- data["snp"]
  rsid <- data["rsid"]
  snpData_alt <- plotCounts(dds = dds, gene = snp,
                            intgroup = c("condition", "donor", "allele"),
                            returnData = TRUE) %>%
    filter(condition == cond) %>%
    group_by(donor) %>% mutate(total = sum(count)) %>% 
    filter(total >= 10) %>%
    filter(all(count >= 2)) %>%
    ungroup() %>%
    filter(allele == "alt") %>%
    mutate(alt_fraction = count/total) %>%
    mutate(snp = snp) %>%
    mutate(rsid = rsid)
  
  # Try a mean
  # mean_altAllele <- data.frame("snp" = snp,
  #                              "rsid" = rsid,
  #                              "mean_alt_Fraction" = mean(snpData_alt$alt_fraction))
  
  
  return(snpData_alt)
  #return(mean_altAllele)
}





## ALL RESULTS
# Read in and format CTL/FNF AI rsids/genes ------------------------------------
CTL <- reformatRSIDdata("data/2023-01-10_AIresCTL_rsids.csv")
FNF <- reformatRSIDdata("data/2023-01-10_AIresFNF_rsids.csv")

# DESeq results --------------------------------------------------------
load("data/2023-01-10_AIresCTL.rda")
load("data/2023-01-10_AIresFNF.rda")
load("data/dds.rda")

# Reformat and join with rsid/gene info
resCTL_df <- reformatDESeqdata(notNA_resCTL, CTL) %>%
  mutate(snp = paste0(.$chr, ":", .$pos, ":", .$ref, ":", .$alt))
resFNF_df <- reformatDESeqdata(notNA_resFNF, FNF) %>%
  mutate(snp = paste0(.$chr, ":", .$pos, ":", .$ref, ":", .$alt))

# Merge and keep distinct SNPs
resALL_df <- bind_rows(resCTL_df, resFNF_df) %>%
  distinct(snp, .keep_all = TRUE)


# Get alt allele fractions for results ------------------------------------

# resALL_altFractions <- apply(resALL_df, 1, get_altFractions, dds = dds) %>%
#   bind_rows() %>%
#   write_csv(paste0("data/", Sys.Date(), "_resALL_altFractions.csv"))



# resCTL_altFractions <- apply(resCTL_df, 1, get_condition_altFractions, dds = dds, cond = "CTL") %>%
#   bind_rows() %>%
#   write_csv(paste0("data/", Sys.Date(), "_resCTL_mean_altFractions.csv"))
# 
# 
# resFNF_altFractions <- apply(resFNF_df, 1, get_condition_altFractions, dds = dds, cond = "FNF") %>%
#   bind_rows() %>%
#   write_csv(paste0("data/", Sys.Date(), "_resFNF_mean_altFractions.csv"))

resCTL_altFractions <- apply(resCTL_df, 1, get_condition_altFractions, dds = dds, cond = "CTL") %>%
  bind_rows() %>%
  write_csv(paste0("data/", Sys.Date(), "_resCTL_altFractions.csv"))


resFNF_altFractions <- apply(resFNF_df, 1, get_condition_altFractions, dds = dds, cond = "FNF") %>%
  bind_rows() %>%
  write_csv(paste0("data/", Sys.Date(), "_resFNF_altFractions.csv"))
  
# Plotting ----------------------------------------------------------------

# ggplot(resALL_altFractions, aes(x = alt_fraction)) +
#   geom_histogram(bins = 40) +
#   scale_x_continuous(limits = c(0,1), expand = c(0, 0), name = "alternative allele fraction",
#                      breaks = c(0, 0.25, 0.50, 0.75, 1),
#                      labels = c(0, 0.25, 0.50, 0.75, 1)) +
#   scale_y_continuous(expand = c(0, 0), name = "Frequency",
#                      limits = c(0, 3000),
#                      breaks = c(0, 500, 1000, 1500, 2000, 2500, 3000),
#                      labels = c(0, 500, 1000, 1500, 2000, 2500, 3000)) +
#   theme_light()
# 
# ggsave(paste0("plots/", Sys.Date(), "_ALLresults_altAlleleFractions_hist.pdf"),
#        units = "in",
#        width = 6.5,
#        height = 5)

ggplot(resCTL_altFractions, aes(x = alt_fraction)) +
  geom_histogram(bins = 40) +
  scale_x_continuous(limits = c(0,1), expand = c(0, 0), name = "alternative allele fraction",
                     breaks = c(0, 0.25, 0.50, 0.75, 1),
                     labels = c(0, 0.25, 0.50, 0.75, 1)) +
  scale_y_continuous(expand = c(0, 0), name = "Frequency",
                     limits = c(0,2000)) +
  theme_light()

ggsave(paste0("plots/", Sys.Date(), "_CTLresults_altAlleleFractions_hist.pdf"),
       units = "in",
       width = 6.5,
       height = 5)
# 
# 
ggplot(resFNF_altFractions, aes(x = alt_fraction)) +
  geom_histogram(bins = 40) +
  scale_x_continuous(limits = c(0,1), expand = c(0, 0), name = "alternative allele fraction",
                     breaks = c(0, 0.25, 0.50, 0.75, 1),
                     labels = c(0, 0.25, 0.50, 0.75, 1)) +
  scale_y_continuous(expand = c(0, 0), name = "Frequency",
                     limits = c(0, 2000)) +
  theme_light()

ggsave(paste0("plots/", Sys.Date(), "_FNFresults_altAlleleFractions_hist.pdf"),
       units = "in",
       width = 6.5,
       height = 5)
# 
# 
# 
# ## SIGNIFICANT RESULTS
# # Read in and format CTL/FNF AI rsids/genes ------------------------------------
# CTL <- reformatRSIDdata("data/2023-01-10_AIsigCTL05_rsids.csv")
# FNF <- reformatRSIDdata("data/2023-01-10_AIsigFNF05_rsids.csv")
# 
# # DESeq results --------------------------------------------------------
# load("data/2023-01-10_AIsigCTL_.05.rda")
# load("data/2023-01-10_AIsigFNF_.05.rda")
# load("data/dds.rda")
# 
# 
# # Format and subset -------------------------------------------------------
# 
# # Reformat and join with rsid/gene info
# sigCTL_df <- reformatDESeqdata(sigCTL_.05, CTL) %>%
#   mutate(snp = paste0(.$chr, ":", .$pos, ":", .$ref, ":", .$alt))
# sigFNF_df <- reformatDESeqdata(sigFNF_.05, FNF) %>%
#   mutate(snp = paste0(.$chr, ":", .$pos, ":", .$ref, ":", .$alt))
# 
# 
# 
# CTL_sigonly <- sigCTL_df[which(!sigCTL_df$rsid %in% sigFNF_df$rsid),]
# FNF_sigonly <- sigFNF_df[which(!sigFNF_df$rsid %in% sigCTL_df$rsid),]
# both_sigonly <- sigCTL_df[which(sigCTL_df$rsid %in% sigFNF_df$rsid),]
# 
# # Get alt allele fractions ------------------------------------------------
# sigCTL_altFractions <- apply(CTL_sigonly, 1, get_condition_altFractions, dds = dds, cond = "CTL") %>%
#   bind_rows() %>%
#   write_csv(paste0("data/", Sys.Date(), "_sigCTL05_altFractions_means.csv"))
# 
# sigFNF_altFractions <- apply(FNF_sigonly, 1, get_condition_altFractions, dds = dds, cond = "FNF") %>%
#   bind_rows() %>%
#   write_csv(paste0("data/", Sys.Date(), "_sigFNF05_altFractions_means.csv"))
# 
# sigboth_altFractions <- apply(both_sigonly, 1, get_altFractions, dds = dds) %>%
#   bind_rows() %>%
#   write_csv(paste0("data/", Sys.Date(), "_sigboth05_altFractions_means.csv"))
# 
# # Plotting ----------------------------------------------------------------
# 
# ggplot(sigCTL_altFractions, aes(x = mean_alt_Fraction)) +
#   geom_histogram(bins = 40) +
#   scale_x_continuous(limits = c(0,1), expand = c(0, 0), name = "alternative allele fraction",
#                      breaks = c(0, 0.25, 0.50, 0.75, 1),
#                      labels = c(0, 0.25, 0.50, 0.75, 1)) +
#   scale_y_continuous(expand = c(0, 0), name = "Frequency",
#                      limits = c(0, 25)) +
#   theme_light()
# 
# ggsave(paste0("plots/", Sys.Date(), "_CTLsigonlyresults_altAlleleFractions_hist.pdf"),
#        units = "in",
#        width = 6.5,
#        height = 5)
# 
# ggplot(sigFNF_altFractions, aes(x = mean_alt_Fraction)) +
#   geom_histogram(bins = 40) +
#   scale_x_continuous(limits = c(0,1), expand = c(0, 0), name = "alternative allele fraction",
#                      breaks = c(0, 0.25, 0.50, 0.75, 1),
#                      labels = c(0, 0.25, 0.50, 0.75, 1)) +
#   scale_y_continuous(expand = c(0, 0), name = "Frequency",
#                      limits = c(0, 15)) +
#   theme_light()
# 
# ggsave(paste0("plots/", Sys.Date(), "_FNFsigonlyresults_altAlleleFractions_hist.pdf"),
#        units = "in",
#        width = 6.5,
#        height = 5)
# 
# 
# ggplot(sigboth_altFractions, aes(x = mean_alt_Fraction)) +
#   geom_histogram(bins = 40) +
#   scale_x_continuous(limits = c(0,1), expand = c(0, 0), name = "alternative allele fraction",
#                      breaks = c(0, 0.25, 0.50, 0.75, 1),
#                      labels = c(0, 0.25, 0.50, 0.75, 1)) +
#   scale_y_continuous(expand = c(0, 0), name = "Frequency") +
#   theme_light()
# 
# ggsave(paste0("plots/", Sys.Date(), "_sigbothresults_altAlleleFractions_hist.pdf"),
#        units = "in",
#        width = 6.5,
#        height = 5)
