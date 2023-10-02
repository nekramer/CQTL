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
  gene_logmu <- design_mat %*% coef_mat[gene,]
  
  geneCounts <- plotCounts(dds, gene = gene,
                           intgroup = "Age",
                           normalized = TRUE,
                           returnData = TRUE) %>%
    rownames_to_column(var = "sample") |> 
    mutate(gene_id = gene,
           logmu = gene_logmu[,1])
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
  gene_clusters <- kmeans(spline_curve_matrix_centered, 2)[["cluster"]]
  
  # Join assigned clusters with gene data
  sig_age_genes$cluster <- gene_clusters
  
  # Get counts and spline fit for genes, split by sample, join with cluster
  sample_gene_countFits <- lapply(sig_age_genes$gene_id,
                                  get_gene_ageCounts, dds_age_lrt) |> 
    bind_rows() |> 
    left_join(sig_age_genes, by = "gene_id") 
  
  
  # Calculate lm per cluster to get general trend (increasing vs decreasing)
  sample_gene_countFits_1 <- sample_gene_countFits |> 
    filter(cluster == 1)
  sample_gene_countFits_2 <- sample_gene_countFits |> 
    filter(cluster == 2)
  cluster_slopes <- data.frame(cluster = c(1, 2),
                               slope = c(coef(lm(count ~ Age, sample_gene_countFits_1))[2],
                                         coef(lm(count ~ Age , sample_gene_countFits_2))[2]))
  # Designate more increasing cluster
  increasingCluster <- cluster_slopes |> 
    filter(slope == max(slope)) |> 
    pull(cluster)
  
  decreasingCluster <- cluster_slopes |> 
    filter(slope == min(slope)) |> 
    pull(cluster)
  
  clusterLabels <- data.frame("cluster" = c(increasingCluster, decreasingCluster),
                              "clusterSign" = c("+", "-"))
  
  # Join with sample gene data
  pval01clusters <- sample_gene_countFits |> 
    left_join(clusterLabels, by = "cluster") |> 
    dplyr::select(-cluster) |> 
    dplyr::rename(cluster = clusterSign) |> 
    write_csv(file = paste0("data/age_de/", filePrefix, "_pval01clusters.csv"))
  
  return(pval01clusters)
  
}

# Wrangle gse -------------------------------------------------------------

# Load gse object
load("data/2023-04-19_gse.rda")

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

ctl_age_pval01clusters <- age_de_analysis(gse, condition = "CTL")

# FNF ---------------------------------------------------------------------

fnf_age_pval01clusters <- age_de_analysis(gse, condition = "FNF")

# Run Homer for GO terms, KEGG pathways, and TF binding motifs ------------

# Compile gene lists for cluster + and cluster -, getting unique genes and
# making sure gene_id is column 6 for run_homer.sh
all_cluster_up <- bind_rows(ctl_age_pval01clusters |> 
                            filter(cluster == "+"),
                          fnf_age_pval01clusters |> 
                            filter(cluster == "+")) |> 
  relocate(gene_id, .after = seqnames) |> 
  distinct(gene_id, .keep_all = TRUE) |> 
  write_csv("data/age_de/pval01_cluster_up.csv")

all_cluster_down <- bind_rows(ctl_age_pval01clusters |> 
                            filter(cluster == "-"),
                          fnf_age_pval01clusters |> 
                            filter(cluster == "-")) |> 
  relocate(gene_id, .after = seqnames) |> 
  distinct(gene_id, .keep_all = TRUE) |> 
  write_csv("data/age_de/pval01_cluster_down.csv")

# Cluster +
system("scripts/run_homer.sh data/age_de/pval01_cluster_up.csv data/homer/homer_age_sig_cluster_up")

# Cluster -
system("scripts/run_homer.sh data/age_de/pval01_cluster_down.csv data/homer/homer_age_sig_cluster_down")



