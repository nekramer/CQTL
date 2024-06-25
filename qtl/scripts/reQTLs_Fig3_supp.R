library(tidyverse)
library(data.table)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(ggvenn)
library(scales)
library(ggpubr)
library(patchwork)
library(googlesheets4)
library(googledrive)
source("../plotting_utils.R")
source("../utils.R")
library(tximeta)
library(org.Hs.eg.db)

# Functions ---------------------------------------------------------------

# Function to convert covariates to factors then numeric and filter out any
# with only 1 unique value to make compatible with Pearson's correlation testing
covariates_to_numeric <- function(cov){
  converted <- cov |> 
    mutate(across(where(is.character), as.factor)) |>
    mutate(across(where(is.Date), as.factor)) |>
    mutate(across(where(is.difftime), as.factor)) |>
    mutate(across(where(is.POSIXt), as.factor)) |> 
    mutate(across(where(is.factor), as.numeric)) |> 
    select_if(~n_distinct(.) > 1)
  
  return(converted)
}

# Function to iterate through covariates and calculate Pearson's correlation and p-value
correlationTests <- function(x, y){
  correlationTestx <- function(x0, y){
    correlationTesty <- function(y0, x0, covarXname){
      covariateName <- colnames(y)[y0]
      result <- cor.test(y[[y0]], x0, method = "pearson")
      return(data.frame("cor" = result$estimate,
                        "pval" = result$p.value,
                        "x" = covarXname,
                        "y" = covariateName))
    }
    # Correlate first df column with each column of the second df
    covariateName <- colnames(x)[x0]
    resY <- sapply(1:ncol(y), correlationTesty, x0 = x[[x0]], 
                   covarXname = covariateName,
                   simplify = FALSE) |> bind_rows()
    
    return(resY)
  }
  
  # Go through each column of the first df
  resX <- sapply(1:ncol(x), correlationTestx, y = y, simplify = FALSE) |> 
    bind_rows() |> 
    remove_rownames()
  
  return(resX)
}

checkVariableCorrelation <- function(var, data){
  
  corVariables <- data |>
    # Subset data where x = var and we're not looking at it's correlation with itself
    filter(x == var) |>
    # Get correlations greater than 0.95
    filter(abs(cor) > 0.95) |>
    pull(y)
  
  return(corVariables)
}

# Define a function to look up the effect size of a lead SNP identified by
# Zeggini in our eqtl datasets
zeggini_effect <- function(eGene_eSNP, pbs_data, fnf_data){
  print(eGene_eSNP)
  eGene_id <- eGene_eSNP[["gene_id"]]
  
  
  snp_chr_lowgrade <- paste0("chr", unlist(str_split(eGene_eSNP[["genotype_id_lowgrade"]], ":"))[1])
  snp_pos_lowgrade <- as.numeric(eGene_eSNP[["hg38pos_lowgrade"]])
  snp_ref_lowgrade <- eGene_eSNP[["REF_lowgrade"]]
  snp_alt_lowgrade <- eGene_eSNP[["ALT_lowgrade"]]
  lowgrade_varIDs <- c(paste0(snp_chr_lowgrade, ":", snp_pos_lowgrade, ":", snp_ref_lowgrade, ":", snp_alt_lowgrade),
                       paste0(snp_chr_lowgrade, ":", snp_pos_lowgrade, ":", snp_alt_lowgrade, ":", snp_ref_lowgrade))
  
  snp_chr_highgrade <- paste0("chr", unlist(str_split(eGene_eSNP[["genotype_id_highgrade"]], ":"))[1])
  snp_pos_highgrade <- as.numeric(eGene_eSNP[["hg38pos_highgrade"]])
  snp_ref_highgrade <- eGene_eSNP[["REF_highgrade"]]
  snp_alt_highgrade <- eGene_eSNP[["ALT_highgrade"]]
  highgrade_varIDs <- c(paste0(snp_chr_highgrade, ":", snp_pos_highgrade, ":", snp_ref_highgrade, ":", snp_alt_highgrade),
                        paste0(snp_chr_highgrade, ":", snp_pos_highgrade, ":", snp_alt_highgrade, ":", snp_ref_highgrade))
  
  if (!is.na(snp_ref_lowgrade)){
    pbs_lowgrade_beta <- pbs_data |> 
      filter(gene_id == eGene_id & 
               variantID %in% lowgrade_varIDs &
               signal == 0) |> 
      # Zeggini reports NES for alt allele, we do minor allele
      mutate(beta = ifelse(ma != snp_alt_lowgrade, -1*beta, beta)) |> 
      pull(beta)
    
    fnf_lowgrade_beta <- fnf_data |> 
      filter(gene_id == eGene_id & variantID %in% lowgrade_varIDs &
               signal == 0) |> 
      # Zeggini reports NES for alt allele, we do minor allele
      mutate(beta = ifelse(ma != snp_alt_lowgrade, -1*beta, beta)) |> 
      pull(beta)
    
    
    if (identical(pbs_lowgrade_beta, logical(0))){
      pbs_lowgrade_beta <- NA
    }
    if (identical(fnf_lowgrade_beta, logical(0))){
      fnf_lowgrade_beta <- NA
    }
  } else {
    pbs_lowgrade_beta <- NA
    fnf_lowgrade_beta <- NA
  }
  
  if (!is.na(snp_ref_highgrade)){
    pbs_highgrade_beta <- pbs_data |> 
      filter(gene_id == eGene_id & variantID %in% highgrade_varIDs &
               signal == 0) |> 
      # Zeggini reports NES for alt allele, we do minor allele
      mutate(beta = ifelse(ma != snp_alt_highgrade, -1*beta, beta)) |> 
      pull(beta)
    
    fnf_highgrade_beta <- fnf_data |> 
      filter(gene_id == eGene_id & variantID %in% highgrade_varIDs &
               signal == 0) |> 
      # Zeggini reports NES for alt allele, we do minor allele
      mutate(beta = ifelse(ma != snp_alt_highgrade, -1*beta, beta)) |> 
      pull(beta)
    
    if (identical(pbs_highgrade_beta, logical(0))){
      pbs_highgrade_beta <- NA
    }
    if (identical(fnf_highgrade_beta, logical(0))){
      fnf_highgrade_beta <- NA
    }
  } else {
    pbs_highgrade_beta <- NA
    fnf_highgrade_beta <- NA
  }
  
  return(tibble(pbs_lowgrade_beta = pbs_lowgrade_beta,
                pbs_highgrade_beta = pbs_highgrade_beta,
                fnf_lowgrade_beta = fnf_lowgrade_beta,
                fnf_highgrade_beta = fnf_highgrade_beta))
  
}

# PEER factor covariate correlation ---------------------------------------

# Read in samplesheets and join for various technical covariates
rnaSamplesheet <- read_csv("data/samplesheet.csv",
                           col_types = "ccccdddTTdddTTccccc") |> 
  dplyr::select(c("Sample", "Donor", "Condition", "Tech_Rep",
                  "Seq_Rep", "RNAextractDate", "RNAextractTime",
                  "RNAextractionKitBatch", "RNAQubit", "RIN",
                  "RNATapeStationDate", "RNAshippedDate", "SequencingBatch")) |> 
  # Convert RNAextractDate and RNAextractTime DateTimes to just dates and times
  mutate(RNAextractDate = as.Date(RNAextractDate),
         RNAextractTime = hms::as_hms(RNAextractTime))

dnaSamplesheet <- read_csv("data/dnaSamplesheet.csv",
                           col_types = c("ccdTddddddddTTcccc")) |>
  dplyr::select(c("Donor", "Date/time DNA extracted", 
                  "DNAReagentBatch", "DNA Qubit conc",
                  "DNA Nanodrop conc.", "A260", "A280", "A270",
                  "A260/A280", "A260/A230", "Date/time passed off to genotyping core",
                  "GenotypingBatch", "Preparer"))


donorSamplesheet <- read_csv("data/donorSamplesheet.csv",
                             col_types = c("ccdcddcDDTddcc")) |> 
  dplyr::select(c("Donor", "Sex", "Age", "OAGradeAvg",
                  "CauseOfDeath", "DeliveryDate", "DateSerumFree",
                  "TimeSerumFree", "FragmentBatch")) |> 
  # Convert TimeSerumFree to just time
  mutate(TimeSerumFree = hms::as_hms(TimeSerumFree))

# Join DNA and Donor samplesheets to RNA by donor
covariates <- rnaSamplesheet |> 
  left_join(dnaSamplesheet, by = "Donor") |> 
  left_join(donorSamplesheet, by = "Donor")

# Read in PEER factors and join with covariate information
ctl_peer_covariate <- read_csv("data/PEERfactors/CTL_PEERfactors_k20.txt") |> 
  left_join(covariates |> filter(Condition == "CTL"), by = "Donor") |> 
  relocate(Sample, Donor)
fnf_peer_covariate <- read_csv("data/PEERfactors/FNF_PEERfactors_k22.txt") |> 
  left_join(covariates |> filter(Condition == "FNF"), by = "Donor") |> 
  relocate(Sample, Donor)

## Pearson's correlation testing

# Convert to numeric
ctl_peer_covariate_numeric <- covariates_to_numeric(ctl_peer_covariate)
fnf_peer_covariate_numeric <- covariates_to_numeric(fnf_peer_covariate)

# Calculate all pairwise correlations

## PBS
ctl_correlations <- correlationTests(ctl_peer_covariate_numeric, 
                                     ctl_peer_covariate_numeric) |> 
  mutate(Condition = "PBS")

# Collapse variables that are highly correlated with each other
ctl_cor_variables <- lapply(unique(ctl_correlations$x), 
                            checkVariableCorrelation, data = ctl_correlations)
names(ctl_cor_variables) <- unique(ctl_correlations$x)
ctl_cor_variables <- ctl_cor_variables |> 
  keep(\(x) length(x) > 1) |> 
  unique()

ctl_correlations_reduced <- ctl_correlations |> 
  filter(!x %in% c("RNAextractDate", "RNAshippedDate", 
                   "SequencingBatch", "RNATapeStationDate") &
           !y %in% c("RNAextractDate", "RNAshippedDate", 
                     "SequencingBatch", "RNATapeStationDate", "Donor", "Sample")) |> 
  mutate(x = ifelse(x == "RNAextractionKitBatch", "RNAbatch", x),
         y = ifelse(y == "RNAextractionKitBatch", "RNAbatch", y)) |> 
  filter(!x %in% c("Date/time DNA extracted", "DNAReagentBatch",
                   "Date/time passed off to genotyping core", "A260", "A280") &
           !y %in% c("Date/time DNA extracted", "DNAReagentBatch",
                     "Date/time passed off to genotyping core", "A260", "A280")) |> 
  mutate(x = ifelse(x == "GenotypingBatch", "DNAbatch", x), 
         y = ifelse(y == "GenotypingBatch", "DNAbatch", y)) |> 
  filter(x != "DateSerumFree" & y != "DateSerumFree") |> 
  mutate(x = ifelse(x == "DeliveryDate", "Donorbatch", x),
         y = ifelse(y == "DeliveryDate", "Donorbatch", y)) |> 
  filter(grepl("PEER", x)) |> 
  filter(!grepl("PEER", y))


## FNF
fnf_correlations <- correlationTests(fnf_peer_covariate_numeric, 
                                     fnf_peer_covariate_numeric) |> 
  mutate(Condition = "FN-f")

# Collapse variables that are highly correlated with each other
fnf_cor_variables <- lapply(unique(fnf_correlations$x), 
                            checkVariableCorrelation, data = fnf_correlations)
names(fnf_cor_variables) <- unique(fnf_correlations$x)
fnf_cor_variables <- fnf_cor_variables |> 
  keep(\(x) length(x) > 1) |> 
  unique()

fnf_correlations_reduced <- fnf_correlations |> 
  filter(!x %in% c("RNAextractDate", "RNAshippedDate", 
                   "SequencingBatch", "RNATapeStationDate") &
           !y %in% c("RNAextractDate", "RNAshippedDate", 
                     "SequencingBatch", "RNATapeStationDate", "Donor", "Sample")) |> 
  mutate(x = ifelse(x == "RNAextractionKitBatch", "RNAbatch", x),
         y = ifelse(y == "RNAextractionKitBatch", "RNAbatch", y)) |> 
  filter(!x %in% c("Date/time DNA extracted", "DNAReagentBatch",
                   "Date/time passed off to genotyping core", "A260", "A280") &
           !y %in% c("Date/time DNA extracted", "DNAReagentBatch",
                     "Date/time passed off to genotyping core", "A260", "A280")) |> 
  mutate(x = ifelse(x == "GenotypingBatch", "DNAbatch", x), 
         y = ifelse(y == "GenotypingBatch", "DNAbatch", y)) |> 
  filter(x != "DateSerumFree" & y != "DateSerumFree") |> 
  mutate(x = ifelse(x == "DeliveryDate", "Donorbatch", x),
         y = ifelse(y == "DeliveryDate", "Donorbatch", y)) |> 
  filter(grepl("PEER", x)) |> 
  filter(!grepl("PEER", y))

# Join for plotting
all_peer_covar_cor <- bind_rows(ctl_correlations_reduced, 
                                fnf_correlations_reduced) |>
  # P-value significance levels
  mutate(pval_sig = case_when(pval < 0.001 ~ "***",
                              pval < 0.01 ~ "**",
                              pval < 0.05 ~ "*")) |> 
  # Condition factor 
  mutate(Condition = factor(Condition, levels = c("PBS", "FN-f"))) |> 
  # Covariate ordering
  mutate(x = factor(x, levels = rev(paste0("PEER", 1:22)))) |> 
  mutate(y = factor(y, levels = c("Sex", "Age",                                  
                                  "CauseOfDeath", "OAGradeAvg", "Donorbatch",
                                  "TimeSerumFree", "FragmentBatch", "RNAbatch",
                                  "RNAQubit", "RNAextractTime", "RIN", 
                                  "Tech_Rep", "Seq_Rep",
                                  "A270", "A260/A280", "A260/A230", 
                                  "DNA Qubit conc", "DNA Nanodrop conc.", 
                                  "DNAbatch", "Preparer")))

peer_covariate_cor_heatmap <- ggplot(all_peer_covar_cor, aes(x = y, y = x, fill = cor)) +
  geom_tile() +
  scale_fill_gradient2(limits = c(-1, 1), low = "#418C82", mid = "white", high = "#C55E2D",
  ) +
  geom_text(aes(label = pval_sig), family = "Helvetica", size = 2, fontface = "bold") +
  facet_wrap(~Condition, scales = "free") +
  guides(fill = guide_colorbar(title = "Pearson's correlation",
                               title.position = "top",
                               title.hjust = 0.5, direction = "horizontal")) +
  theme(axis.text.x = element_text(angle = 45, size = 5, vjust = 1, hjust = 1,
                                   color = "black"),
        axis.text.y = element_text(size = 6, color = "black"),
        panel.background = element_blank(),
        strip.background = element_blank(),
        text = element_text(family = "Helvetica"),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        legend.title = element_text(size = 6, color = "black"),
        legend.text = element_text(size = 5, color = "black"),
        legend.position = "bottom",
        legend.margin = margin(t = -20),
        strip.text = element_text(size = 10))

save(peer_covariate_cor_heatmap, 
     file = "plots/reQTLs_Fig3_supp/peer_covariate_cor_heatmap.rda")



# PEER/covariate eGene analysis -------------------------------------------

resultPath <- "data/eqtl/qtl_supp/"
CTL_resultFiles <- list.files(resultPath, 
                              pattern = paste0("^CTL_PEER_k.*_perm1Mb_sig\\.csv$"))

FNF_resultFiles <- list.files(resultPath, 
                              pattern = paste0("^FNF_PEER_k.*_perm1Mb_sig\\.csv$"))

# Get number of significant eGenes in PBS and FN-f for each PEER/covariate dataset
peer_sig_eGenes <- list()

# Add final CTL sig results
Nk <- 20
batchGroup <- ""
num_eGenes <- read_csv("data/eqtl/CTL_PEER_k20_genoPC_perm1Mb_sig_rsID.csv") |> 
  nrow()

peer_sig_eGenes[["CTL_final"]] <- data.frame("Condition" = "PBS",
                                             "PEER" = Nk,
                                             "eGenes" = num_eGenes,
                                             "batchGroup" = batchGroup)

for (file in CTL_resultFiles){
  
  Nk <- as.numeric(gsub("k", "", unlist(str_split(file, "_"))[3]))
  RNAbatch <- "RNAKitBatch" %in% unlist(str_split(file, "_"))
  DNAbatch <- "DNAKitBatch" %in% unlist(str_split(file, "_"))
  sequencingBatch <- "RNASequencingBatch" %in% unlist(str_split(file, "_"))
  genotypingBatch <- "genoBatch" %in% unlist(str_split(file, "_"))
  
  batchGroup <- paste(c("RNAKitBatch", 
                        "DNAKitBatch", 
                        "RNAsequencingBatch", 
                        "genoBatch")[c(RNAbatch,
                                       DNAbatch,
                                       sequencingBatch,
                                       genotypingBatch)], collapse = "_")
  
  num_eGenes <- read_csv(paste0(resultPath, file)) |> 
    nrow()
  
  
  peer_sig_eGenes[[file]] <- data.frame("Condition" = "PBS",
                                        "PEER" = Nk,
                                        "eGenes" = num_eGenes,
                                        "batchGroup" = batchGroup)
}

# Add final FNF sig results
Nk <- 22
batchGroup <- ""
num_eGenes <- read_csv("data/eqtl/FNF_PEER_k22_genoPC_perm1Mb_sig_rsID.csv") |> 
  nrow()

peer_sig_eGenes[["FNF_final"]] <- data.frame("Condition" = "FN-f",
                                             "PEER" = Nk,
                                             "eGenes" = num_eGenes,
                                             "batchGroup" = batchGroup)

for (file in FNF_resultFiles){
  
  Nk <- as.numeric(gsub("k", "", unlist(str_split(file, "_"))[3]))
  RNAbatch <- "RNAKitBatch" %in% unlist(str_split(file, "_"))
  DNAbatch <- "DNAKitBatch" %in% unlist(str_split(file, "_"))
  sequencingBatch <- "RNASequencingBatch" %in% unlist(str_split(file, "_"))
  genotypingBatch <- "genoBatch" %in% unlist(str_split(file, "_"))
  
  batchGroup <- paste(c("RNAKitBatch", 
                        "DNAKitBatch", 
                        "RNAsequencingBatch", 
                        "genoBatch")[c(RNAbatch,
                                       DNAbatch,
                                       sequencingBatch,
                                       genotypingBatch)], collapse = "_")
  
  num_eGenes <- read_csv(paste0(resultPath, file)) |>
    nrow()
  
  peer_sig_eGenes[[file]] <- data.frame("Condition" = "FN-f",
                                        "PEER" = Nk,
                                        "eGenes" = num_eGenes,
                                        "batchGroup" = batchGroup)
}

peer_sig_eGenes <- peer_sig_eGenes |> 
  bind_rows() |> 
  mutate(batchGroup = ifelse(batchGroup == "", "No batches", batchGroup)) |> 
  mutate(batchGroup = ifelse(batchGroup == "RNAKitBatch_DNAKitBatch_RNAsequencingBatch_genoBatch", 
                             "RNAKitBatch_DNAKitBatch\nRNAsequencingBatch_genoBatch", batchGroup)) |> 
  mutate(batchGroup = factor(batchGroup, 
                             levels = c("DNAKitBatch",
                                        "RNAKitBatch",
                                        "RNAKitBatch_DNAKitBatch",
                                        "RNAKitBatch_DNAKitBatch\nRNAsequencingBatch_genoBatch",
                                        "No batches"))) |> 
  mutate(Condition = factor(Condition, levels = c("PBS", "FN-f")))


batchLabels_all <- peer_sig_eGenes |> 
  filter(PEER == 50) |> 
  filter(batchGroup != "No batches") |> 
  mutate(y = case_when(batchGroup == "RNAKitBatch" ~ eGenes + 10,
                       batchGroup == "RNAKitBatch_DNAKitBatch" ~ eGenes - 100,
                       batchGroup == "DNAKitBatch" ~ eGenes - 50,
                       batchGroup == "RNAKitBatch_DNAKitBatch\nRNAsequencingBatch_genoBatch" ~ eGenes))

batchLabels_none <- peer_sig_eGenes |> 
  filter(PEER == 50) |> 
  filter(batchGroup == "No batches") |> 
  mutate(y = eGenes + 50)

peer_sig_egenes_lines <- ggplot(data = peer_sig_eGenes, 
       aes(x = PEER, y = eGenes, color = batchGroup)) +
  geom_vline(data = filter(peer_sig_eGenes, Condition == "PBS"),
             aes(xintercept = 20), color = "grey", lty = 2) +
  geom_text(data = filter(peer_sig_eGenes, Condition == "PBS"),
            aes(x = 21, y = 2800, label = "20"), color = "grey", hjust = 0,
            family = "Helvetica") +
  geom_vline(data = filter(peer_sig_eGenes, Condition == "FN-f"),
             aes(xintercept = 22), color = "grey", lty = 2) +
  geom_text(data = filter(peer_sig_eGenes, Condition == "FN-f"),
            aes(x = 23, y = 2800, label = "22"), color = "grey", hjust = 0,
            family = "Helvetica") +
  geom_line(lwd = 0.75) +
  coord_cartesian(clip = "off") +
  geom_text(data = batchLabels_all, aes(x = Inf, y = y, label = batchGroup),
            family = "Helvetica", size = 3, hjust = 0, vjust = 1) +
  geom_text(data = batchLabels_none, aes(x = Inf, y = y, label = batchGroup),
            family = "Helvetica", size = 3, fontface = "bold", hjust = 0, vjust = 1) +
  facet_wrap(~Condition, scales = "free") +
  scale_y_continuous(limits = c(0, 3000)) +
  scale_color_manual(values = c("grey75", "grey75", "grey75", "grey75", "#136079")) +
  xlab("Number of PEER factors") +
  ylab("Number of eGenes") +
  theme_custom_scatterplot() +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 10),
        axis.line = element_line(color = "black"),
        axis.title = element_text(color = "black", size = 8),
        axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        plot.margin = unit(c(1, 10, 1, 1), "lines"),
        legend.position = "none",
        panel.spacing = unit(10, "lines"))

save(peer_sig_egenes_lines, file = "plots/reQTLs_Fig3_supp/peer_sig_egenes_lines.rda")

# All PBS/FNF eQTL venn diagram -------------------------------------------

CTL_eQTLs <- list()
FNF_eQTLs <- list()

for (chrom in 1:22){
  ctl_chrom_qtls <- fread(paste0("data/eqtl/qtl_nom/CTL_PEER_k20_genoPC_allSignals_nom1Mb_MAFs_chr",
                                 chrom, ".csv"), data.table = FALSE) |>
    filter(nom_sig == 1)
  CTL_eQTLs[[chrom]] <- ctl_chrom_qtls

  fnf_chrom_qtls <- fread(paste0("data/eqtl/qtl_nom/FNF_PEER_k22_genoPC_allSignals_nom1Mb_MAFs_chr",
                                 chrom, ".csv"), data.table = FALSE) |>
    filter(nom_sig == 1)
  FNF_eQTLs[[chrom]] <- fnf_chrom_qtls

}

CTL_eQTLs <- bind_rows(CTL_eQTLs) |>
  mutate(eGene_variant = paste0(gene_id, "_", variantID))
FNF_eQTLs <- bind_rows(FNF_eQTLs) |>
  mutate(eGene_variant = paste0(gene_id, "_", variantID))

pbs_fnf_eqtls <- tibble(values = unique(c(CTL_eQTLs$eGene_variant,
                                          FNF_eQTLs$eGene_variant))) |>
  mutate(PBS = values %in% CTL_eQTLs$eGene_variant,
         FNF = values %in% FNF_eQTLs$eGene_variant)

pbs_fnf_all_eqtls_venn <- ggplot(pbs_fnf_eqtls, aes(A = PBS, B = FNF)) +
  geom_venn(set_names = c("PBS", "FN-f"),
            fill_color = c(log2fcColors[["-"]], log2fcColors[["+"]]),
            stroke_color = NA, auto_scale = TRUE, show_percentage = FALSE,
            text_size = 6, set_name_size = 6) +
  coord_fixed()  +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank())

pbs_fnf_all_eqtls_venn  <- venn_font(pbs_fnf_all_eqtls_venn , font = "Helvetica")

save(pbs_fnf_all_eqtls_venn, 
     file = "plots/reQTLs_Fig3_supp/pbs_fnf_all_eqtls_venn.rda")

# eQTL supp figure --------------------------------------------------------
pdf(file = "plots/reQTLs_Fig3_supp/SupFig3.pdf",
    width = 11, height = 8.75)
pageCreate(width = 11, height = 8.75, showGuides = FALSE)

### A - PEER correlation with covariates

plotGG(peer_covariate_cor_heatmap, x = 0.1, y = 0, width = 5.5, height = 4.5)
plotText("A", x = 0.1, y = 0.1, just = c("left", "top"), fontfamily = "Helvetica",
         fontsize = 14, fontface = "bold")

## B - eGene selection plots
plotText("B", x = 0.1, y = 4.5, just = c("left", "top"), fontfamily = "Helvetica",
         fontsize = 14, fontface = "bold")
plotGG(peer_sig_egenes_lines, x = 0, y = 4.25, width = 11, height =4.55)

## C - venn diagram of all eGene-eSNP pairs
plotGG(pbs_fnf_all_eqtls_venn, x = 6.1, y = 0, width = 4.5, height = 4.5)
plotText("C", x = 5.75, y = 0.1, just = c("left", "top"), fontfamily = "Helvetica",
         fontsize = 14, fontface = "bold")
dev.off()

# Venn diagram and effect eGene overlap with Steinberg et al. 2021 -------------------------

# Combine results from lowgrade and highgrade cartilage results
lowgrade_Zeggini <- read_csv("data/Steinberg_2021/processed/eQTL_LowGradeCartilage_perm_sig_lead_hg38.csv",
                             col_types = "cdcddddddddcddddddddddddddccc") |> 
  separate_wider_delim(cols = "phenotype_id", delim = "_", 
                       names = c("gene_symbol", "gene_id"))
highgrade_Zeggini <- read_csv("data/Steinberg_2021/processed/eQTL_HighGradeCartilage_perm_sig_lead_hg38.csv",
                              col_types = "cdcddddddddcddddddddddddddccc") |> 
  separate_wider_delim(cols = "phenotype_id", delim = "_", 
                       names = c("gene_symbol", "gene_id"))
lowgrade_highgrade_Zeggini_eGenes <- unique(c(lowgrade_Zeggini$gene_id,
                                              highgrade_Zeggini$gene_id))

# PBS and FN-f eGenes
PBS_eGenes <- read_csv("data/eqtl/CTL_PEER_k20_genoPC_cond1Mb_topSignals_rsID.csv")
FNF_eGenes <- read_csv("data/eqtl/FNF_PEER_k22_genoPC_cond1Mb_topSignals_rsID.csv")
PBS_FNF_eGenes <- unique(c(PBS_eGenes$gene_id,
                           FNF_eGenes$gene_id))

# Create tibble of groups
pbsfnf_highlowgrade_eGenes <- tibble(values = unique(c(PBS_FNF_eGenes,
                                                       lowgrade_highgrade_Zeggini_eGenes))) |> 
  mutate(Kramer = values %in% PBS_FNF_eGenes,
         Zeggini = values %in% lowgrade_highgrade_Zeggini_eGenes)

kramer_zeggini_venn <- ggplot(pbsfnf_highlowgrade_eGenes , 
                              aes(A = Zeggini, B = Kramer)) +
  geom_venn(set_names = c("Steinberg et al.", "Kramer et al."),
            fill_color = c("#4A8D9D", "#A9D48A"),
            stroke_color = NA, auto_scale = TRUE, show_percentage = FALSE, 
            text_size = 3, set_name_size = 3.5) +
  coord_fixed() +
  theme(panel.background = element_rect(fill = 'transparent'),
        plot.background = element_rect(fill = 'transparent'),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank())
kramer_zeggini_venn <- venn_font(kramer_zeggini_venn, font = "Helvetica")

save(kramer_zeggini_venn, file = "plots/reQTLs_Fig3_supp/kramer_zeggini_venn.rda")


## Format table of shared and not shared eGenes/leads
kramer_pbs <- read_csv("data/eqtl/CTL_PEER_k20_genoPC_cond1Mb_topSignals_rsID.csv", 
         col_select = c("gene_id", "variantID", "signal")) |> 
  filter(signal == 0) |> 
  dplyr::select(-signal)
kramer_fnf <- read_csv("data/eqtl/FNF_PEER_k22_genoPC_cond1Mb_topSignals_rsID.csv", 
                       col_select = c("gene_id", "variantID", "signal")) |> 
  filter(signal == 0) |> 
  dplyr::select(-signal)

kramer_pbs_fnf <- full_join(kramer_pbs, kramer_fnf, by = c("gene_id"),
                            suffix = c("_Kramer_PBS", "_Kramer_FNF"))

zeggini_lowgrade <- read_csv("data/Steinberg_2021/processed/eQTL_LowGradeCartilage_perm_sig_lead_hg38.csv",
                             col_types = "cdcddddddddcddddddddddddddccc") |> 
  separate_wider_delim(cols = "phenotype_id", delim = "_", 
                       names = c("gene_symbol", "gene_id")) |> 
  separate_wider_delim(cols = "genotype_id", delim = ":",
                       names = c("chr", "hg19pos")) |> 
  mutate(variantID = paste0("chr", chr, ":", hg38pos, ":", REF, ":", ALT)) |> 
  dplyr::select(gene_id, variantID) |> 
  distinct()

zeggini_highgrade <- read_csv("data/Steinberg_2021/processed/eQTL_HighGradeCartilage_perm_sig_lead_hg38.csv",
                              col_types = "cdcddddddddcddddddddddddddccc") |> 
  separate_wider_delim(cols = "phenotype_id", delim = "_", 
                       names = c("gene_symbol", "gene_id")) |> 
  separate_wider_delim(cols = "genotype_id", delim = ":",
                       names = c("chr", "hg19pos")) |> 
  mutate(variantID = paste0("chr", chr, ":", hg38pos, ":", REF, ":", ALT)) |> 
  dplyr::select(gene_id, variantID) |> 
  distinct()

zeggini_lowgrade_highgrade <- full_join(zeggini_lowgrade, zeggini_highgrade,
                                        by = c("gene_id"),
                                        suffix = c("_Steinberg_lowgrade", "_Steinberg_highgrade"))

mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
kramer_zeggini_all <- full_join(kramer_pbs_fnf, zeggini_lowgrade_highgrade,
                                by = c("gene_id")) |> 
  left_join(getBM(filters = "ensembl_gene_id", attributes = c("ensembl_gene_id","hgnc_symbol"),
                  values = kramer_zeggini_all$gene_id, mart = mart),
            by = join_by("gene_id" == "ensembl_gene_id")) |> 
  dplyr::rename(gene_symbol = hgnc_symbol) |> 
  relocate(gene_symbol, .after = "gene_id") |> 
  rowwise() |> 
  mutate(study = 
           case_when(any(!is.na(variantID_Kramer_PBS), 
                         !is.na(variantID_Kramer_FNF)) &
                       any(!is.na(variantID_Steinberg_lowgrade),
                           !is.na(variantID_Steinberg_highgrade)) ~ "both",
                     any(!is.na(variantID_Kramer_PBS),
                         !is.na(variantID_Kramer_FNF)) &
                       all(is.na(variantID_Steinberg_lowgrade),
                           is.na(variantID_Steinberg_highgrade)) ~ "Kramer",
                     all(is.na(variantID_Kramer_PBS),
                         is.na(variantID_Kramer_FNF)) &
                       any(!is.na(variantID_Steinberg_lowgrade),
                           !is.na(variantID_Steinberg_highgrade)) ~ "Steinberg")) |> 
  arrange(variantID_Kramer_PBS)


write_csv(kramer_zeggini_all, file = "tables/SupTable8.csv")


#### Direction of effect of shared study eGenes
shared_kramer_zeggini <- pbsfnf_highlowgrade_eGenes |> 
  filter(Kramer == TRUE & Zeggini == TRUE) |> 
  pull(values)

shared_kramer_egenes <- full_join(PBS_eGenes |> 
            filter(gene_id %in% shared_kramer_zeggini) |> 
            dplyr::select(gene_id) |> 
            mutate(PBS = "PBS"),
          FNF_eGenes |> 
            filter(gene_id %in% shared_kramer_zeggini) |> 
            dplyr::select(gene_id) |> 
            mutate(FNF = "FN-f"), by = "gene_id")


shared_zeggini_egenes <- full_join(lowgrade_Zeggini |> 
                                     filter(gene_id %in% shared_kramer_zeggini) |> 
                                     dplyr::select(gene_id, gene_symbol, PERM_slope, genotype_id, hg38pos, 
                                                   REF, ALT),
                                   highgrade_Zeggini |> 
                                     filter(gene_id %in% shared_kramer_zeggini) |> 
                                     dplyr::select(gene_id, gene_symbol, PERM_slope, genotype_id, hg38pos, 
                                                   REF, ALT),
                                   by = "gene_id", suffix = c("_lowgrade", "_highgrade"),
                                   relationship = "many-to-many") |> 
  # Join with our study conditions
  left_join(shared_kramer_egenes, by = "gene_id") |> 
  distinct()

## Read in nominal data for all chromosomes in both conditions, subsetting for
# shared zeggini egene variants
pbs_nom_data <- list()
fnf_nom_data <- list()
for (chr in 1:22){
  pbs_chr <- fread(paste0("data/eqtl/qtl_nom/CTL_PEER_k20_genoPC_allSignals_nom1Mb_MAFs_chr", 
               chr, ".csv"),
        data.table = FALSE) |> 
    filter(gene_id %in% shared_zeggini_egenes$gene_id)
  
  fnf_chr <- fread(paste0("data/eqtl/qtl_nom/FNF_PEER_k22_genoPC_allSignals_nom1Mb_MAFs_chr", 
                          chr, ".csv"),
                   data.table = FALSE) |> 
    filter(gene_id %in% shared_zeggini_egenes$gene_id)
  
  pbs_nom_data[[paste0("chr", chr)]] <- pbs_chr
  fnf_nom_data[[paste0("chr", chr)]] <- fnf_chr
}

pbs_nom_data_zeggini <- pbs_nom_data |> bind_rows()
fnf_nom_data_zeggini <- fnf_nom_data |> bind_rows()
## Look up SNP effects in our study

shared_zeggini_egene_effects <- apply(shared_zeggini_egenes, 1, zeggini_effect, pbs_data = pbs_nom_data_zeggini, 
                                      fnf_data = fnf_nom_data_zeggini)
  
save(shared_zeggini_egene_effects, file = "data/Steinberg_2021/processed/shared_zeggini_egene_effects.rda")

shared_zeggini_egene_effects <- bind_rows(shared_zeggini_egene_effects)


shared_zeggini_egenes_kramer_effects <- bind_cols(shared_zeggini_egenes,
                                                  shared_zeggini_egene_effects) |> 
  mutate(dir_lowgrade = ifelse(PERM_slope_lowgrade < 0, "down", "up"),
         dir_highgrade = ifelse(PERM_slope_highgrade < 0, "down", "up"),
         dir_pbs_lowgrade = ifelse(pbs_lowgrade_beta < 0, "down", "up"),
         dir_pbs_highgrade = ifelse(pbs_highgrade_beta < 0, "down", "up"),
         dir_fnf_lowgrade = ifelse(fnf_lowgrade_beta < 0, "down", "up"),
         dir_fnf_highgrade = ifelse(fnf_highgrade_beta < 0, "down", "up"))


# Lowgrade comparisons to PBS and FN-f
lowgrade_pbs <- ggplot(shared_zeggini_egenes_kramer_effects,
       mapping = aes(x = pbs_lowgrade_beta, y = PERM_slope_lowgrade)) +
  geom_vline(xintercept = 0, color = "grey20", linetype = 2) +
  geom_hline(yintercept = 0, color = "grey20", linetype = 2) +
  geom_point(color = "grey20", size = 0.5) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  stat_cor(label.y = 2, family = "Helvetica", size = 3) +
  scale_x_continuous(name = "beta (Kramer et al.; PBS)") +
  scale_y_continuous(name = "beta (Steinberg et al., 2021; low-grade cartilage)") +
  theme_custom_scatterplot() +
  theme(axis.title = element_text(family = "Helvetica",
                                  size = 8),
        axis.text = element_text(family = "Helvetica", size = 8))

lowgrade_fnf <- ggplot(shared_zeggini_egenes_kramer_effects,
       mapping = aes(x = fnf_lowgrade_beta, y = PERM_slope_lowgrade)) +
  geom_vline(xintercept = 0, color = "grey20", linetype = 2) +
  geom_hline(yintercept = 0, color = "grey20", linetype = 2) +
  geom_point(color = "grey20", size = 0.5) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  stat_cor(label.y = 2, family = "Helvetica", size = 3) +
  scale_x_continuous(name = "beta (Kramer et al.; FN-f)") +
  scale_y_continuous(name = "beta (Steinberg et al., 2021; low-grade cartilage)") +
  theme_custom_scatterplot() +
  theme(axis.title = element_text(family = "Helvetica",
                                  size = 8),
        axis.text = element_text(family = "Helvetica", size = 8))


# Highgrade comparisons to PBS and FN-f
highgrade_pbs <- ggplot(shared_zeggini_egenes_kramer_effects,
       mapping = aes(x = pbs_highgrade_beta, y = PERM_slope_highgrade)) +
  geom_vline(xintercept = 0, color = "grey20", linetype = 2) +
  geom_hline(yintercept = 0, color = "grey20", linetype = 2) +
  geom_point(color = "grey20", size = 0.5) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  stat_cor(label.y = 2, family = "Helvetica", size = 3) +
  scale_x_continuous(name = "beta (Kramer et al.; PBS)") +
  scale_y_continuous(name = "beta (Steinberg et al., 2021; high-grade cartilage)") +
  theme_custom_scatterplot() +
  theme(axis.title = element_text(family = "Helvetica",
                                  size = 8),
        axis.text = element_text(family = "Helvetica", size = 8))

highgrade_fnf <- ggplot(shared_zeggini_egenes_kramer_effects,
       mapping = aes(x = fnf_highgrade_beta, y = PERM_slope_highgrade)) +
  geom_vline(xintercept = 0, color = "grey20", linetype = 2) +
  geom_hline(yintercept = 0, color = "grey20", linetype = 2) +
  geom_point(color = "grey20", size = 0.5) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  stat_cor(label.y = 2, family = "Helvetica", size = 3) +
  scale_x_continuous(name = "beta (Kramer et al.; FN-f)") +
  scale_y_continuous(name = "beta (Steinberg et al., 2021; high-grade cartilage)") +
  theme_custom_scatterplot() +
  theme(axis.title = element_text(family = "Helvetica",
                                  size = 8),
        axis.text = element_text(family = "Helvetica", size = 8))

kramer_steinberg_overlap_scatterplots <- 
  (lowgrade_pbs + lowgrade_fnf) / (highgrade_pbs + highgrade_fnf) + 
  plot_annotation(theme = theme(panel.background = 
                                  element_rect(fill = "transparent", 
                                               color = "transparent"),
                                plot.background = 
                                  element_rect(fill = "transparent", 
                                               color = "transparent")))
save(kramer_steinberg_overlap_scatterplots, file = "plots/reQTLs_Fig3_supp/kramer_steinberg_overlap_scatterplots.rda")


# Steinberg/Zeggini overlap supp figure 3 ---------------------------------------------------

pdf(file = "plots/reQTLs_Fig3_supp/SupFig5.pdf",
    width = 8.5, height = 5.5)
pageCreate(width = 8.5, height = 5.5, showGuides = FALSE)

### A - Venn diagram of eGene overlap
plotText("A", x = 0.1, y = 0.1, just = c("left", "top"), fontfamily = "Helvetica",
         fontsize = 14, fontface = "bold")
plotGG(kramer_zeggini_venn, x = -0.3, y = -0.2, width = 3.5, height = 3.5)

### B - Scatterplots of betas for steinberg lead snp of shared genes
plotText("B", x = 2.8, y = 0.1, just = c("left", "top"), fontfamily = "Helvetica",
         fontsize = 14, fontface = "bold")
plotGG(kramer_steinberg_overlap_scatterplots, x = 3, y = 0.1, width = 5.5, height = 5.5)

dev.off()





# Format eQTL and condition-specific QTL tables --------------------------

load("data/rna/qtl_gse.rda")
geneInfo <- rowRanges(gse) |> 
  as_tibble() |> 
  rename(gene_end = end,
         gene_start = start)

### Standard eQTLs
## PBS

#Perm pass
PBS_perm <- read_csv("data/eqtl/CTL_PEER_k20_genoPC_perm1Mb_sig_rsID.csv",
                     col_select = c("gene_id", "qval")) |> 
  rename(eGene_qval = qval)

# Conditional pass
PBS_eQTLs <- read_csv("data/eqtl/CTL_PEER_k20_genoPC_cond1Mb_topSignals_rsID.csv") |> 
  dplyr::select(gene_id, gene_symbol, gene_chr, gene_start,
                gene_strand, rsID, variantID, variant_start, signal) |> 
  dplyr::rename(gene_name = gene_symbol,
                gene_tss = gene_start,
                variant_pos = variant_start) |>
  left_join(PBS_perm, by = "gene_id") |> 
  left_join(geneInfo |> 
              dplyr::select(gene_id, gene_end), by = "gene_id") |> 
  relocate(gene_end, .after = gene_tss)

PBS_nom_data <- list()
for (chr in 1:22){
  chrom_nomData <- fread(paste0("data/eqtl/qtl_nom/CTL_PEER_k20_genoPC_allSignals_nom1Mb_MAFs_chr",
                            chr,".csv"), data.table = FALSE, 
                         select = c("gene_id", "variantID", "num_cis_variants",
                                    "beta", "beta_se", "pval_nominal_threshold",
                                    "nom_pval",
                                    "ma", "maf", "signal")) |> 
    filter(gene_id %in% PBS_eQTLs$gene_id) |> 
    rename(minor_allele = ma,
           MAF = maf,
           eGene_nominal_threshold = pval_nominal_threshold)
  PBS_nom_data[[paste0("chr", chr)]] <- chrom_nomData
}
PBS_nom_data <- bind_rows(PBS_nom_data)


PBS_eQTLs_final <- PBS_eQTLs |> 
  left_join(PBS_nom_data, by = c("gene_id", "variantID", "signal")) |> 
  relocate(num_cis_variants, .after = gene_strand) |> 
  relocate(eGene_nominal_threshold, .after = MAF) |> 
  relocate(eGene_qval, .after = eGene_nominal_threshold) |> 
  relocate(signal, .after = eGene_qval) |> 
  group_by(signal) |> 
  arrange(signal) |> 
  arrange(gene_name, .by_group = TRUE) |> 
  ungroup()

write_csv(PBS_eQTLs_final, file = "tables/SupTable7_PBS.csv")

## FNF

# Perm pass
FNF_perm <- read_csv("data/eqtl/FNF_PEER_k22_genoPC_perm1Mb_sig_rsID.csv",
                     col_select = c("gene_id", "qval")) |> 
  rename(eGene_qval = qval)

# Conditional pass
FNF_eQTLs <- read_csv("data/eqtl/FNF_PEER_k22_genoPC_cond1Mb_topSignals_rsID.csv") |> 
  dplyr::select(gene_id, gene_symbol, gene_chr, gene_start,
                gene_strand, rsID, variantID, variant_start, signal) |> 
  dplyr::rename(gene_name = gene_symbol,
                gene_tss = gene_start,
                variant_pos = variant_start) |>
  left_join(FNF_perm, by = "gene_id") |> 
  left_join(geneInfo |> 
              dplyr::select(gene_id, gene_end), by = "gene_id") |> 
  relocate(gene_end, .after = gene_tss)

FNF_nom_data <- list()
for (chr in 1:22){
  chrom_nomData <- fread(paste0("data/eqtl/qtl_nom/FNF_PEER_k22_genoPC_allSignals_nom1Mb_MAFs_chr",
                                chr,".csv"), data.table = FALSE, 
                         select = c("gene_id", "variantID", "num_cis_variants",
                                    "beta", "beta_se", "pval_nominal_threshold",
                                    "nom_pval",
                                    "ma", "maf", "signal")) |> 
    filter(gene_id %in% FNF_eQTLs$gene_id) |> 
    rename(minor_allele = ma,
           MAF = maf,
           eGene_nominal_threshold = pval_nominal_threshold)
  FNF_nom_data[[paste0("chr", chr)]] <- chrom_nomData
}
FNF_nom_data <- bind_rows(FNF_nom_data)

FNF_eQTLs_final <- FNF_eQTLs |> 
  left_join(FNF_nom_data, by = c("gene_id", "variantID", "signal")) |> 
  relocate(num_cis_variants, .after = gene_strand) |> 
  relocate(eGene_nominal_threshold, .after = MAF) |> 
  relocate(eGene_qval, .after = eGene_nominal_threshold) |> 
  relocate(signal, .after = eGene_qval) |> 
  group_by(signal) |> 
  arrange(signal) |> 
  arrange(gene_name, .by_group = TRUE) |> 
  ungroup()

write_csv(FNF_eQTLs_final, file = "tables/SupTable7_FNF.csv")

ss <- gs4_create(name = "SupTable7")
write_sheet(PBS_eQTLs_final,
            ss, sheet = "PBS")
write_sheet(FNF_eQTLs_final,
            ss, sheet = "FN-f")

drive_mv(file = "SupTable7", path = as_dribble("CQTL paper/Figures and Tables"))

### Condition-specific eQTLs
## PBS
PBS_cond_highconf <- read_csv("data/reqtl/CTL_sig01_beta_donor_reQTLs_PEER_k20_genoPC.csv",
                              col_select = "gene_id")

PBS_cond <- read_csv("data/reqtl/CTL_sig01_reQTLs_PEER_k20_genoPC.csv") |> 
  dplyr::select(gene_id, interaction_pval) |> 
  left_join(PBS_eQTLs_final |> 
              filter(signal == 0), by = "gene_id") |> 
  mutate(high_conf = ifelse(gene_id %in% PBS_cond_highconf$gene_id, "high confidence", NA)) |> 
  dplyr::select(-signal,-num_cis_variants,-eGene_nominal_threshold,-eGene_qval, -nom_pval) |> 
  relocate(interaction_pval, .after = MAF) |> 
  rename(PBS_beta = beta,
         PBS_beta_se = beta_se)

pbs_cond_fnf_beta <- list()
for (chr in 1:22){
  chrom_nomData <- fread(paste0("data/eqtl/qtl_nom/FNF_PEER_k22_genoPC_allSignals_nom1Mb_MAFs_chr",
                                chr,".csv"), data.table = FALSE, 
                         select = c("gene_id", "variantID",
                                    "beta", "beta_se", "signal")) |> 
    filter(gene_id %in% PBS_cond$gene_id & signal == 0) |> 
    dplyr::select(-signal) |> 
    rename(FNF_beta = beta,
           FNF_beta_se = beta_se)
  pbs_cond_fnf_beta[[paste0("chr", chr)]] <- chrom_nomData
  
}
pbs_cond_fnf_beta <- bind_rows(pbs_cond_fnf_beta)

PBS_cond_final <- PBS_cond |> 
  left_join(pbs_cond_fnf_beta, by = c("gene_id", "variantID")) |> 
  relocate(FNF_beta, .after = PBS_beta_se) |> 
  relocate(FNF_beta_se, .after = FNF_beta) |> 
  group_by(high_conf) |> 
  arrange(high_conf) |> 
  arrange(gene_name, .by_group = TRUE) |> 
  ungroup()

write_csv(PBS_cond_final, file = "tables/SupTable9_PBS.csv")

## FNF
FNF_cond_highconf <- read_csv("data/reqtl/FNF_sig01_beta_donor_reQTLs_PEER_k22_genoPC.csv",
                              col_select = "gene_id")

FNF_cond <- read_csv("data/reqtl/FNF_sig01_reQTLs_PEER_k22_genoPC.csv") |> 
  dplyr::select(gene_id, interaction_pval) |> 
  left_join(FNF_eQTLs_final |> 
              filter(signal == 0), by = "gene_id") |> 
  mutate(high_conf = ifelse(gene_id %in% FNF_cond_highconf$gene_id, "high confidence", NA)) |> 
  dplyr::select(-signal,-num_cis_variants,-eGene_nominal_threshold,-eGene_qval, -nom_pval) |> 
  relocate(interaction_pval, .after = MAF) |> 
  rename(FNF_beta = beta,
         FNF_beta_se = beta_se)

fnf_cond_pbs_beta <- list()
for (chr in 1:22){
  chrom_nomData <- fread(paste0("data/eqtl/qtl_nom/CTL_PEER_k20_genoPC_allSignals_nom1Mb_MAFs_chr",
                                chr,".csv"), data.table = FALSE, 
                         select = c("gene_id", "variantID",
                                    "beta", "beta_se", "signal")) |> 
    filter(gene_id %in% FNF_cond$gene_id & signal == 0) |> 
    dplyr::select(-signal) |> 
    rename(PBS_beta = beta,
           PBS_beta_se = beta_se)
  fnf_cond_pbs_beta[[paste0("chr", chr)]] <- chrom_nomData
  
}
fnf_cond_pbs_beta <- bind_rows(fnf_cond_pbs_beta)

FNF_cond_final <- FNF_cond |> 
  left_join(fnf_cond_pbs_beta, by = c("gene_id", "variantID")) |> 
  relocate(PBS_beta, .before = FNF_beta) |> 
  relocate(PBS_beta_se, .before = FNF_beta) |> 
  group_by(high_conf) |> 
  arrange(high_conf) |> 
  arrange(gene_name, .by_group = TRUE) |> 
  ungroup()

write_csv(FNF_cond_final, file = "tables/SupTable9_FNF.csv")

ss <- gs4_create(name = "SupTable9")
write_sheet(PBS_cond_final,
            ss, sheet = "PBS")
write_sheet(FNF_cond_final,
            ss, sheet = "FN-f")

drive_mv(file = "SupTable9", path = as_dribble("CQTL paper/Figures and Tables"))

