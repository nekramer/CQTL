library(tidyverse)
library(corrplot)
library(limma)
# Functions ---------------------------------------------------------------


pca_to_plot <- function(condition){
  samplesheet <- read_csv("../samplesheet.csv") %>%
    # Get distinct samples
    distinct(Sample, .keep_all = TRUE)
  
  CPMadjTMM_invNorm <- read_csv(paste0("../output/normquant/", condition, "_CPMadjTMM_invNorm.bed"))
  CPMadjTMM_invNorm_info <- CPMadjTMM_invNorm %>% 
    dplyr::select(-seqnames, -start, -end, -strand, -gene_name) %>%
    t() %>% as.data.frame() %>% row_to_names(row_number = 1) %>%
    rownames_to_column("Sample") %>%
    left_join(samplesheet) 
  # Calculate variances to remove columns where variance is 0
  variances <- CPMadjTMM_invNorm %>% 
    dplyr::select(-gene_id, -seqnames, -start, -end, -strand, -gene_name) %>% 
    t() %>% apply(2, var, na.rm = TRUE)
  if (length(which(variances == 0) > 0)){
    
    CPMadjTMM_invNorm_pca <- CPMadjTMM_invNorm %>% 
      dplyr::select(-gene_id, -seqnames, -start, -end, -strand, -gene_name) %>% 
      t() %>% as.data.frame() %>% dplyr::select(-which(variances == 0)) %>%
      prcomp(center = TRUE, scale = TRUE)
    
  } else {
    CPMadjTMM_invNorm_pca <- CPMadjTMM_invNorm %>% 
      dplyr::select(-gene_id, -seqnames, -start, -end, -strand, -gene_name) %>% 
      t() %>%
      prcomp(center = TRUE, scale = TRUE)
  }
  
  CPMadjTMM_invNorm_info$PC1 <- CPMadjTMM_invNorm_pca$x[,1]
  CPMadjTMM_invNorm_info$PC2 <- CPMadjTMM_invNorm_pca$x[,2]
  assign(paste0(condition, "_CPMadjTMM_invNorm_info"), CPMadjTMM_invNorm_info, envir = globalenv())
  assign(paste0(condition, "_PCs"), CPMadjTMM_invNorm_pca$x[,1:10], envir = globalenv())
}

corInfo <- function(condition){
  
  data <- get(paste0(condition,"_CPMadjTMM_invNorm_info"), envir = globalenv())
  info <- data[,c("RNAQubit.x", "Sex", "Age", "Race", "OAGradeAvg", "CauseOfDeath", "Fragment Batch", "Donor.y", "Time.y", "Tech_Rep.y", "Seq_Rep.y", "RNAextractionKitBatch", "RNAQubit.y", "RIN.y")] %>% 
    mutate(across(where(is_character),as_factor))  %>% 
    mutate(across(where(is.factor),as.numeric))
  info$RIN.y <- as.numeric(unlist(info$RIN.y))
  assign(paste0(condition, "_info"), info, envir = globalenv())
}

correlationTests_p <- function(x, y){
  
  correlationTestx <- function(x0, y){
    
    correlationTesty <- function(x0, y0){
      result <- cor.test(x0, y0, method = "pearson")
      return(result$p.value)
    }
    
    apply(y, 2, correlationTesty, y0 = x0)
    
  }
  
  apply(x, 2, correlationTestx, y = y)
  
  
}

correlationTests_cor <- function(x, y){
  
  correlationTestx <- function(x0, y){
    
    correlationTesty <- function(x0, y0){
      result <- cor.test(x0, y0, method = "pearson")
      return(result$estimate)
    }
    
    apply(y, 2, correlationTesty, y0 = x0)
    
  }
  
  apply(x, 2, correlationTestx, y = y)
  
}

# Spreadsheet information -------------------------------------------------

donors <- read_sheet(ss = "https://docs.google.com/spreadsheets/d/1JwLw9D6rMqhHC9BPrZebAN40Wojo-CqbMdiIPXAzkLo/edit#gid=1699779981",
                     sheet = "Donors")
donors$OAGradeAvg <- as.numeric(unlist(donors$OAGradeAvg))
donors$`Fragment Batch` <- as.numeric(unlist(donors$`Fragment Batch`))

rna <- read_sheet(ss = "https://docs.google.com/spreadsheets/d/1JwLw9D6rMqhHC9BPrZebAN40Wojo-CqbMdiIPXAzkLo/edit#gid=1699779981",
                  sheet = "RNAExtractionsLibraries") %>% filter(NYGCQC == "PASS") %>% filter(!is.na(Sequencing_Directory)) #%>% 
                #select(Sample, Donor, Condition, Time, Tech_Rep, Seq_Rep, RNAQubit, RIN, RNAextractionKitBatch)
rna$RNAextractionKitBatch <- as.numeric(unlist(rna$RNAextractionKitBatch))
rna$Tech_Rep <- as.numeric(unlist(rna$Tech_Rep))
rna$Seq_Rep <- as.numeric(unlist(rna$Seq_Rep))

# No batch normalization --------------------------------------------------

lapply(c("CTL", "FNF"), pca_to_plot)


# Tried coloring by Sex, Race, RNAExtractionKitBatch, Age and showed no clear patterns
CTL_CPMadjTMM_invNorm_info <- CTL_CPMadjTMM_invNorm_info %>% select(-RNAextractionKitBatch) %>% left_join(donors, by = "Donor") %>% left_join(rna, by = c("Sample")) %>% # Get distinct samples
  distinct(Sample, .keep_all = TRUE)

FNF_CPMadjTMM_invNorm_info <- FNF_CPMadjTMM_invNorm_info %>% select(-RNAextractionKitBatch) %>% left_join(donors, by = "Donor") %>% left_join(rna, by = c("Sample")) %>% # Get distinct samples
  distinct(Sample, .keep_all = TRUE)


# CTL plot
ggplot(data = CTL_CPMadjTMM_invNorm_info, mapping = aes(x = PC1, y = PC2, color = Donor)) +
  geom_point() + theme_minimal()


# FNF plot
ggplot(data = FNF_CPMadjTMM_invNorm_info, mapping = aes(x = PC1, y = PC2, color = Donor)) +
  geom_point() + theme_minimal()


# Correlation testing PC's and technical variables ----------------------------

lapply(c("CTL", "FNF"), corInfo)


# Testing correlation of first 10 PC's with various technical variables
CTL_cor <- correlationTests_cor(CTL_PCs, CTL_info) %>% na.omit()
CTL_cor_p <- correlationTests_p(CTL_PCs, CTL_info) %>% na.omit()

FNF_cor <- correlationTests_cor(FNF_PCs, FNF_info) %>% na.omit()
FNF_cor_p <- correlationTests_p(FNF_PCs, FNF_info) %>% na.omit()

corrplot(CTL_cor, p.mat = CTL_cor_p, sig.level = 0.05)
corrplot(FNF_cor, p.mat = FNF_cor_p, sig.level = 0.05)


# Pearson's correlation between donor replicates -------------------------------

rep_samplesheet <- read_csv("../samplesheet3.csv") %>%
  # Get distinct samples
  distinct(Sample, .keep_all = TRUE)

repadjTMM_invNorm <- read_csv("../output/normquant/ALLwithreps_CPMadjTMM_invNorm.bed") %>%
  # Grab rep donors
  select(matches("AM7205|AM7204|AM7224")) %>%
  relocate(matches("FNF"), .after = matches("CTL")) %>%
  relocate(matches("AM7204"), .before = matches("AM7205|AM7224")) %>%
  relocate(matches("AM7205"), .before =  matches("AM7224")) %>%
  relocate(ends_with("_1"), .before = ends_with("_2"))

Repcorrelations <- cor(repadjTMM_invNorm)

corrplot(Repcorrelations)


Repcorrelationsdf <- as.data.frame(Repcorrelations)
controlCorrelations <- Repcorrelationsdf %>% 
  select(contains("CTL")) %>% 
  rownames_to_column(var = "Sample") %>%
  filter(grepl("CTL", Sample) == TRUE)

CTL_cor_vector <- c()
for (i in 1:nrow(controlCorrelations)){
  # Get sample
  sample <- controlCorrelations[i,"Sample"]
  # Only use first reps to find second reps
  if (endsWith(sample, "_1")){
    # Get donor 
    donor <- str_split(sample, "_", simplify = TRUE)[1,2]
    
    # Determine what rep sample would be
    repSample <- paste0("CQTL_", donor, "_R_CTL_18_2")
    
    # Grab correlation value from this column
    correlation <- controlCorrelations[i, repSample]

    CTL_cor_vector[i] <- correlation
    
  }
}

fnfCorrelations <- Repcorrelationsdf %>% 
  select(contains("FNF")) %>% 
  rownames_to_column(var = "Sample") %>%
  filter(grepl("FNF", Sample) == TRUE)

FNF_cor_vector <- c()
for (i in 1:nrow(fnfCorrelations)){
  # Get sample
  sample <- fnfCorrelations[i,"Sample"]
  # Only use first reps to find second reps
  if (endsWith(sample, "_1")){
    # Get donor 
    donor <- str_split(sample, "_", simplify = TRUE)[1,2]
    
    # Determine what rep sample would be
    repSample <- paste0("CQTL_", donor, "_R_FNF_18_2")
    
    # Grab correlation value from this column
    correlation <- fnfCorrelations[i, repSample]
    
    FNF_cor_vector[i] <- correlation
    
  }
}

# Combine data for ggplot boxplot
CTL_cor <- data.frame("cor" = CTL_cor_vector[which(!is.na(CTL_cor_vector))], "condition" = "CTL")
FNF_cor <- data.frame("cor" = FNF_cor_vector, "condition" = "FNF")

all_cors <- bind_rows(CTL_cor, FNF_cor)

ggplot(data = all_cors, mapping = aes(x = condition, y = cor)) +geom_boxplot() + theme_minimal()
