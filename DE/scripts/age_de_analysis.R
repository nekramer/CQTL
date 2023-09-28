library(DESeq2)
library(tidyverse)
library(plyranges)
library(splines)


# Functions ---------------------------------------------------------------

# Function to extract vector of spline fit values for a gene
get_gene_splineFit <- function(gene, dds){
  # Get coef matrix and design matrix for fitted spline
  coef_mat <- coef(dds)
  design_mat <- model.matrix(design(dds), colData(dds))
  
  # Spline fit row
  splineFit <- as.data.frame(t(design_mat %*% coef_mat[gene,]))
  return(splineFit) 
}

# Function to get normalized sample counts and spline fit values for a gene in 
# a dataframe
get_gene_ageCounts <- function(gene, dds){
  # Get coef matrix and design matrix for fitted spline
  coef_mat <- coef(dds)
  design_mat <- model.matrix(design(dds), colData(dds))
  geneCounts <- plotCounts(dds, gene = gene,
                           intgroup = "Age",
                           normalized = TRUE,
                           returnData = TRUE) %>%
    rownames_to_column(var = "sample") |> 
    mutate(gene_id = gene,
           logmu = design_mat %*% coef_mat[gene,])
  return(geneCounts)
}

# Function to carry out age de analysis for a condition
age_de_analysis <- function(gse, condition){
  
  if (condition == "CTL"){
    gse_analysis <- gse[, gse$Condition == "CTL"]
    dds_name <- "dds_age_ctl_lrt"
    filePrefix <- "ctl_age"
  } else if (condition == "FNF"){
    gse_analysis <- gse[, gse$Condition == "FNF"]
    dds_name <- "dds_age_fnf_lrt"
    filePrefix <- "fnf_age"
  }
  
  # Build DESeq object
  dds_age <- DESeqDataSet(gse_analysis, 
                          design = ~Sex + Race + ns(Age, df = 5))
  
  # Filter out lowly expressed genes
  keep <- rowSums(counts(dds_age) >= 10) >= ceiling(nrow(colData(gse_analysis))*0.5)
  dds_age <- dds_age[keep,]
  
  # Perform LRT with reduced model without Age
  dds_age_lrt <- DESeq(dds_age, test = "LRT", reduced = ~Sex + Race)
  
  # Rename and save dds
  assign(dds_name, dds_age_lrt)
  save(list = dds_name, file = paste0("data/age_de/", dds_name, ".rda"))
  
  # Get significant genes and format
  sig_age_genes <- results(dds_age_lrt, format = "GRanges") |> 
    filter(padj < 0.01) |> 
    plyranges::names_to_column("gene_id") |> 
    as.data.frame() |> 
    inner_join(as.data.frame(rowData(gse_analysis)) |> 
                 dplyr::select(c("gene_id", "symbol", "tx_ids")),
               by = "gene_id") |> 
    makeGRangesFromDataFrame(keep.extra.columns = TRUE) |> 
    keepStandardChromosomes(pruning.mode = "coarse") |> 
    as.data.frame()
  
  # Build spline curve matrix of significant hits
  spline_curve_matrix <- lapply(sig_age_genes$gene_id, 
                                get_gene_splineFit, 
                                dds_age_lrt) |> bind_rows() 
  rownames(spline_curve_matrix) <- sig_age_genes$gene_id
  
  # Center matrix
  spline_curve_matrix_centered <- 
    spline_curve_matrix - rowMeans(spline_curve_matrix)
  
  # Perform clustering
  gene_clusters <- kmeans(spline_curve_matrix, 2)[["cluster"]]
  
  # Join clusters with gene data
  sig_age_genes$cluster <- gene_clusters 
  
  # Get counts and spline fit for genes, split by sample, and join with gene 
  # cluster information
  lapply(sig_age_genes$gene_id,
         get_gene_age_Counts, dds_age_lrt) |> 
    bind_rows() |> 
    left_join(sig_age_genes, by = "gene_id") |> 
    write_csv(file = paste0("data/age_de/", filePrefix, "_pval01clusters.csv"))
  
}

# Wrangle gse -------------------------------------------------------------

# Load gse object
load("data/2023-04-19_gse.rda")
#load("data/2023-09-27_gse.rda")

# Read in donorSamplesheet for additional donor info
donorSamplesheet <- read_csv("data/donorSamplesheet.csv") |> 
  mutate(Race = replace_na(Race, "Unknown"))

# Join gse colData with donorSamplesheet
colData(gse) <- as(left_join(as.data.frame(colData(gse)),
                             donorSamplesheet[,c("Donor", "Sex", "Age", "Race")],
                             by = "Donor"), "DataFrame")

# Add colnames based on sample name/Age
colnames(gse) <- paste0(colData(gse)[,"names"], "_", colData(gse)[,"Age"])

# CTL ---------------------------------------------------------------------

age_de_analysis(gse, condition = "CTL")

# FNF ---------------------------------------------------------------------

age_de_analysis(gse, condition = "FNF")



