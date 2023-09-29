library(DESeq2)
library(tidyverse)
library(plyranges)

# Functions ---------------------------------------------------------------

sex_de_analysis <- function(gse, condition){
  
  if (condition == "CTL"){
    gse_analysis <- gse[, gse$Condition == "CTL"]
    dds_name <- "dds_sex_ctl"
    filePrefix <- "ctl_sex"
  } else if (condition == "FNF"){
    gse_analysis <- gse[, gse$Condition == "FNF"]
    dds_name <- "dds_sex_fnf"
    filePrefix <- "fnf_sex"
  }
  
  # Build DESeq object
  dds_sex <- DESeqDataSet(gse_analysis, design = ~Race + Age_group + Sex)
  
  # Filter out lowly expressed genes
  keep <- rowSums(counts(dds_sex) >= 10) >= ceiling(nrow(colData(gse_analysis))*0.5)
  dds_sex <- dds_sex[keep,]
  
  # Fit model
  dds_sex <- DESeq(dds_sex)
  
  # Rename and save dds
  assign(dds_name, dds_sex)
  save(list = dds_name, file = paste0("data/sex_de/", dds_name, ".rda"))
 
  # Shrink l2fc
  sex_shrink <- lfcShrink(dds_sex,
                          coef = "Sex_M_vs_F", format = "GRanges") |>
    plyranges::names_to_column("gene_id")
  
  # Join results with gene info
  sex_shrink <-
    inner_join(x = as.data.frame(sex_shrink),
               y = as.data.frame(rowData(gse_analysis)) %>%
                 dplyr::select(c("gene_id", "symbol", "tx_ids")),
               by = "gene_id") %>%
    makeGRangesFromDataFrame(keep.extra.columns = TRUE) %>%
    keepStandardChromosomes(pruning.mode = "coarse") %>%
    as.data.frame()
   
  
  ## Save l2fc-shrunken results
  write_csv(sex_shrink, file = paste0("data/sex_de/", filePrefix, "_shrink.csv"))
  
  # Get significant genes and write to file
  sex_shrink |> 
    filter(padj < 0.01) |> 
    write_csv(paste0("data/sex_de/", filePrefix, "DE_pval01.csv"))
}

# gse wrangling -----------------------------------------------------------

# Load gse object
load("data/2023-04-19_gse.rda")

# Read in donorSamplesheet for additional donor info
donorSamplesheet <- read_csv("data/donorSamplesheet.csv") |> 
  mutate(Race = replace_na(Race, "Unknown"))

# Join gse colData with donorSamplesheet
colData(gse) <- as(left_join(as.data.frame(colData(gse)),
                             donorSamplesheet[,c("Donor", "Sex", "Age", "Race")],
                             by = "Donor"), "DataFrame")

## Put Age into categories
gse$Age_group <- cut(gse$Age, breaks = seq(30, 90, 10),
                     labels = c("31-40", "41-50", "51-60", "61-70", "71-80", "81-90"))

# Convert colData to factors
colData(gse)[] <- lapply(colData(gse), factor)

# Add colnames based on sample name/Sex
colnames(gse) <- paste0(colData(gse)[,"names"], "_", colData(gse)[,"Sex"])

# CTL ---------------------------------------------------------------------

sex_de_analysis(gse, condition = "CTL")

# FNF ---------------------------------------------------------------------

sex_de_analysis(gse, condition = "FNF")

# RESPONSE TO TREATMENT ---------------------------------------------------

# Create ind.n in colData to number within each group
colData(gse) <- as(colData(gse) |>
                     as.data.frame() |> 
                     group_by(Sex) |> 
                     mutate(ind.n = as.factor(dense_rank(Donor))), "DataFrame")


# Create design matrix and remove columns for unbalanced M/F interactions
designMatrix <- model.matrix(~Sex + Sex:ind.n + Sex:Condition, colData(gse))
designMatrix <- designMatrix[,-(which(colSums(designMatrix) == 0))]

# Build DESeq object
dds_sex_treatment_response <- DESeqDataSet(gse, 
                                           design = designMatrix)

# Filter out lowly expressed genes
keep <- rowSums(counts(dds_sex_treatment_response) >= 10) >= ceiling(nrow(colData(gse))*0.5)
dds_sex_treatment_response <- dds_sex_treatment_response[keep,]

# Fit model
dds_sex_treatment_response <- DESeq(dds_sex_treatment_response)
## Save object
save(dds_sex_treatment_response, file = "data/sex_de/dds_sex_treatment_response.rda")

# Get results and write to file
sexDE_treatmenteffect_pval01 <- results(dds_sex_treatment_response,
                                        contrast = list("SexF.ConditionFNF",
                                                        "SexM.ConditionFNF"), 
                                        format = "GRanges") |> 
  plyranges::names_to_column("gene_id") |>
  as.data.frame() |> 
  inner_join(y = as.data.frame(rowData(gse)) |> 
               dplyr::select(c("gene_id", "symbol", "tx_ids")),
             by = "gene_id") |> 
  makeGRangesFromDataFrame(keep.extra.columns = TRUE) |> 
  keepStandardChromosomes(pruning.mode = "coarse") |> 
  as.data.frame() |> 
  filter(padj < 0.01) |> 
  write_csv(file = "data/sex_de/sexDE_treatmenteffect_pval01.csv")

