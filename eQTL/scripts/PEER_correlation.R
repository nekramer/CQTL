library(tidyverse)

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

# Covariate parsing -------------------------------------------------------

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
  dplyr::select(c("Donor", "Sex", "Age", "Race", "OAGradeAvg",
                  "CauseOfDeath", "DeliveryDate", "DateSerumFree",
                  "TimeSerumFree", "FragmentBatch")) |> 
  # Convert TimeSerumFree to just time
  mutate(TimeSerumFree = hms::as_hms(TimeSerumFree))

# Join DNA and Donor samplesheets to RNA by donor
covariates <- rnaSamplesheet |> 
  left_join(dnaSamplesheet, by = "Donor") |> 
  left_join(donorSamplesheet, by = "Donor")

# Read in PEER factors and join with covariate information
ctl_peer_covariate <- read_csv("data/CTL_PEERfactors_k20.txt") |> 
  left_join(covariates |> filter(Condition == "CTL"), by = "Donor") |> 
  relocate(Sample, Donor)
fnf_peer_covariate <- read_csv("data/FNF_PEERfactors_k20.txt") |> 
  left_join(covariates |> filter(Condition == "FNF"), by = "Donor") |> 
  relocate(Sample, Donor)


# Pearson's correlation testing -------------------------------------------

# Convert to numeric
ctl_peer_covariate_numeric <- covariates_to_numeric(ctl_peer_covariate)
fnf_peer_covariate_numeric <- covariates_to_numeric(fnf_peer_covariate)

# Calculate all pairwise correlations

## CTL
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
                     "SequencingBatch", "RNATapeStationDate")) |> 
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
         y = ifelse(y == "DeliveryDate", "Donorbatch", y))
  
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
                     "SequencingBatch", "RNATapeStationDate")) |> 
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
         y = ifelse(y == "DeliveryDate", "Donorbatch", y))


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
  mutate(x = factor(x, levels = c("Sample", "Donor", "Sex", "Age", "Race",
                                  "CauseOfDeath", "OAGradeAvg", "Donorbatch",
                                  "TimeSerumFree", "FragmentBatch", "RNAbatch",
                                  "RNAQubit", "RIN", 
                                  "Tech_Rep", "Seq_Rep",
                                  "A270", "A260/A280", "A260/A230", 
                                  "DNA Qubit conc", "DNA Nanodrop conc.", 
                                  "DNAbatch", "Preparer",
                                  paste0("PEER", 1:20)))) |> 
  mutate(y = factor(y, levels = c("Sample", "Donor", "Sex", "Age", "Race",
                                  "CauseOfDeath", "OAGradeAvg", "Donorbatch",
                                  "TimeSerumFree", "FragmentBatch", "RNAbatch",
                                  "RNAQubit", "RIN", 
                                  "Tech_Rep", "Seq_Rep",
                                  "A270", "A260/A280", "A260/A230", 
                                  "DNA Qubit conc", "DNA Nanodrop conc.", 
                                  "DNAbatch", "Preparer",
                                  paste0("PEER", 1:20))))


# Filter for lower triangular data
mask <- with(all_peer_covar_cor, as.integer(x) >= as.integer(y))
all_peer_covar_cor_lt <- all_peer_covar_cor[which(mask), ]

ggplot(all_peer_covar_cor_lt, aes(x = x, y = y, fill = cor)) +
  geom_tile() +
  scale_fill_gradient2(limits = c(-1, 1), low = "#418C82", mid = "white", high = "#C55E2D",
                      ) +
  geom_text(aes(label = pval_sig), family = "Helvetica", size = 2, fontface = "bold") +
  facet_wrap(~Condition) +
  coord_equal() +
  guides(fill = guide_colorbar(title = "Pearson's correlation",
                               title.position = "top",
                               title.hjust = 0.5, direction = "horizontal")) +
  theme(axis.text.x = element_text(angle = 45, size = 6, vjust = 1, hjust = 1,
                                   color = "black"),
        axis.text.y = element_text(size = 6, color = "black"),
        panel.background = element_blank(),
        strip.background = element_blank(),
        text = element_text(family = "Helvetica"),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        legend.title = element_text(size = 6, color = "black"),
        legend.text = element_text(size = 5, color = "black"),
        legend.position = c(0.1,0.7),
        strip.text = element_text(size = 12))

ggsave(filename = "plots/PEERfactors/PEER_covariate_correlation.pdf",
       width = 16, height = 10)
