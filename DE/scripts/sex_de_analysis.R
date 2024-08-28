library(DESeq2)
library(tidyverse)
library(plyranges)
source("../utils.R")

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
  dds_sex <- DESeqDataSet(gse_analysis, design = ~Ancestry + Age_group + Sex)
  
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
load("data/2023-10-03_gse.rda")
# Remove geno contaminated donor
gse <- gse[, gse$Donor != "AM7352"]

donorSamplesheet <- read_csv("data/donorSamplesheet.csv") |>
  dplyr::select(Donor, Sex, Age) |>
  # Read in and join ancestries determined through genotyping pca
  left_join(read_csv("data/CQTL_COA_01_GDA8_COA2_01_COA3_01_GDA8_COA4_COA5_COA6_COA7_predictedAncestry.csv") |>
              separate_wider_delim(cols = "Donor", delim = "_",
                                   names = c(NA, "Donor", NA), too_many = "drop"), by = "Donor") |>
  dplyr::rename(Ancestry = Predicted_Ancestry)

# Join gse colData with donorSamplesheet
colData(gse) <- as(left_join(as.data.frame(colData(gse)),
                             donorSamplesheet[,c("Donor", "Sex", "Age", "Ancestry")],
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

# GO/KEGG on autosomal genes ----------------------------------------------

ctl_sig_genes <- read_csv("/proj/phanstiel_lab/Data/processed/CQTL/rna/CQTL_AM7180_AM7634/sex_de/ctl_sexDE_pval01.csv") |> 
  dplyr::select(seqnames, gene_id)
fnf_sig_genes <- read_csv("/proj/phanstiel_lab/Data/processed/CQTL/rna/CQTL_AM7180_AM7634/sex_de/fnf_sexDE_pval01.csv") |> 
  dplyr::select(seqnames, gene_id)
# Get union list

bind_rows(ctl_sig_genes, fnf_sig_genes) |> 
  distinct() |> 
  # Filter for autososmes
  filter(!seqnames %in% c("X", "Y")) |> 
  dplyr::select(-seqnames) |> 
  write_delim(file = "data/sex/autosome_sex_de_genes.txt",
              col_names = FALSE)

## GO
sex_autosome_go <- 
  read_delim("data/sex/autosome_sex_homer/biological_process.txt") |> 
  mutate(pval = exp(1)^logP) |> 
  filter(pval < 0.01)
sex_autosome_go_red <- reduceGO(sex_autosome_go,
                               category = "sex_autosome")

sex_autosome_go_red |> 
  dplyr::select(-`Entrez Gene IDs`, -pval, -logP, -size, -termUniqueness, -termUniquenessWithinCluster,
                -termDispensability, -category) |> 
  relocate(`-log10pval`, .after = Enrichment) |> 
  arrange(desc(`-log10pval`)) |> 
  write_csv(file = "sex_autosome_GO_sigres.csv")
## KEGG

sex_autosome_kegg <- read_delim("data/sex/autosome_sex_homer/kegg.txt") |> 
  mutate(pval = exp(1)^logP) |> 
  filter(pval < 0.01) |> 
  distinct(Term, .keep_all = TRUE) |> 
  mutate(`-log10pval` = -log10(pval))



# variability across sex de genes -----------------------------------------

ctl_sig_genes <- read_csv("/proj/phanstiel_lab/Data/processed/CQTL/rna/CQTL_AM7180_AM7634/sex_de/ctl_sexDE_pval01.csv")
fnf_sig_genes <- read_csv("/proj/phanstiel_lab/Data/processed/CQTL/rna/CQTL_AM7180_AM7634/sex_de/fnf_sexDE_pval01.csv")

union_sig_genes <- union(ctl_sig_genes$gene_id, fnf_sig_genes$gene_id)

load("/proj/phanstiel_lab/Data/processed/CQTL/rna/CQTL_AM7180_AM7634/sex_de/dds_sex_ctl.rda")
load("/proj/phanstiel_lab/Data/processed/CQTL/rna/CQTL_AM7180_AM7634/sex_de/dds_sex_fnf.rda")

sex_gene_counts <- list()
for (geneid in union_sig_genes){
  
  ctl_gene_counts <- get_gene_sex_Counts(gene = geneid, dds = dds_sex_ctl) |> 
    mutate(condition = "PBS")
  fnf_gene_counts <- get_gene_sex_Counts(gene = geneid, dds = dds_sex_fnf) |> 
    mutate(condition = "FN-f")
  all_gene_counts <- bind_rows(ctl_gene_counts, fnf_gene_counts)
  sex_gene_counts[[geneid]] <- all_gene_counts
}

sex_gene_counts <- bind_rows(sex_gene_counts)



sex_intersect <- intersect(ctl_sig_genes$gene_id,
                           fnf_sig_genes$gene_id)

ctl_only <- ctl_sig_genes$gene_id[!ctl_sig_genes$gene_id %in% sex_intersect]
fnf_only <- fnf_sig_genes$gene_id[!fnf_sig_genes$gene_id %in% sex_intersect]


sex_gene_counts |> 
  mutate(gene_category = case_when(gene_id %in% sex_intersect ~ "PBS and FN-f",
                                   gene_id %in% ctl_only ~ "PBS only",
                                   gene_id %in% fnf_only ~ "FN-f only"),
         condition = factor(condition, levels = c("PBS", "FN-f"))) |> 
  filter(!is.na(gene_category)) |> 
  ggplot(aes(x = gene_id, y = log2(count),
                            color = condition)) +
  geom_point(position = position_dodge(width = 0.7), size = 0.25) +
  geom_boxplot(fill = NA, outlier.shape = NA) +
  facet_wrap(vars(gene_category), nrow = 3, ncol = 1, drop = TRUE, scales = "free") +
  scale_color_manual(values = c(log2fcColors[["-"]], log2fcColors[["+"]])) +
  theme_minimal() +
  theme(axis.text.x = element_blank())


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



