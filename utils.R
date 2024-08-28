# Function to reduce similar GO terms
reduceGO <- function(goterms, category, ont = "BP", threshold = 0.7){
  
  # Calculate similarity matrix for GO terms based on ontology
  simMatrix <- calculateSimMatrix(goterms$TermID,
                                  orgdb = "org.Hs.eg.db",
                                  ont = ont,
                                  method = "Rel")
  
  # Create named vector of scores
  # Higher is better -> -log10 transforming p-values
  scores <- setNames(-log10(goterms$pval), goterms$TermID)
  
  # Group GO terms based on similarity threshold
  reducedTerms <- reduceSimMatrix(simMatrix,
                                  scores,
                                  threshold = threshold,
                                  orgdb = "org.Hs.eg.db") |> 
    dplyr::rename(`-log10pval` = score) 
  
  # Join grouped terms with original data
  joined_reducedTerms <- left_join(goterms, reducedTerms, 
                                   by = join_by("TermID" == "go", "Term" == "term")) |>
    dplyr::select(-cluster) |> 
    relocate(parent, .after = Term) |> 
    relocate(parentTerm, .after = parent) |> 
    dplyr::rename(parentTermID = parent) |> 
    mutate("category" = category)
  
  return(joined_reducedTerms)
}

# Function to get normalized donor gene counts with Condition
get_gene_condition_Counts <- function(gene, dds){
  geneCounts <- DESeq2::plotCounts(dds, gene = gene,
                                   intgroup = "Condition",
                                   normalized = TRUE,
                                   returnData = TRUE) |> 
    remove_rownames() |> 
    mutate(gene_id = gene)
  
  return(geneCounts)
}


# Function to get normalized donor gene counts with sex
get_gene_sex_Counts <- function(gene, dds){
  
  geneCounts <- tryCatch({DESeq2::plotCounts(dds, gene = gene, 
                                   intgroup = "Sex", 
                                   normalized = TRUE,
                                   returnData = TRUE) |> 
    remove_rownames() |> 
    mutate(gene_id = gene)}, error = function(cond) {tibble(row_names = rownames(colData(dds)),
                                                            count = rep(NA, nrow(colData(dds))),
                                                            Sex = colData(dds)[,"Sex"]) |> 
        column_to_rownames(var = "row_names")})
  
  return(geneCounts)
}

# Function to get the normalized donor gene counts with age
get_gene_age_Counts <- function(gene, dds){
  
  
  geneCounts <- tryCatch({DESeq2::plotCounts(dds, gene = gene, 
                                   intgroup = "Age", 
                                   normalized = TRUE,
                                   returnData = TRUE) %>%
    remove_rownames() %>%
    mutate(gene_id = gene)}, error = function(cond) {tibble(row_names = rownames(colData(dds)),
                                                            count = rep(NA, nrow(colData(dds))),
                                                            Age = colData(dds)[,"Age"]) |> 
        column_to_rownames(var = "row_names")})
  
  return(geneCounts)
}

# Function to obtain the minor allele for a given variant ID
get_minor_allele <- function(varID, data){
  # Get chrom
  chrom <- unlist(str_split(varID, ":"))[1]
  
  #eGene
  eGene <- data |> 
    filter(variantID == varID) |> 
    pull(gene_id) |> 
    unique()
  
  # signal number
  signal_no <- data |> 
    filter(variantID == varID) |> 
    pull(signal) |> 
    unique()
  
  # Read in nominal data for that variant's chromosome
  minor_allele <- fread(paste0("data/eqtl/qtl_nom/CTL_PEER_k20_genoPC_allSignals_nom1Mb_MAFs_", chrom , ".csv"), 
                        data.table = FALSE) |> 
    filter(variantID == varID & gene_id == eGene & signal == signal_no) |> 
    pull(ma)
  
  return(minor_allele)
}

# Function to get nominal p-values and betas for a given variant ID
get_pvals_betas <- function(varID, data){
  # Get chrom
  chrom <- unlist(str_split(varID, ":"))[1]
  
  #eGene
  eGene <- data |> 
    filter(variantID == varID) |> 
    pull(gene_id) |> 
    unique()
  
  # signal number
  signal_no <- data |> 
    filter(variantID == varID) |> 
    pull(signal) |> 
    unique()
  
  pbs_data <- fread(paste0("data/eqtl/qtl_nom/CTL_PEER_k20_genoPC_allSignals_nom1Mb_MAFs_", chrom , ".csv"), 
                    data.table = FALSE) |> 
    filter(variantID == varID & gene_id == eGene & signal == signal_no)
  
  fnf_data <- fread(paste0("data/eqtl/qtl_nom/FNF_PEER_k22_genoPC_allSignals_nom1Mb_MAFs_", chrom , ".csv"), 
                    data.table = FALSE) |> 
    filter(variantID == varID & gene_id == eGene & signal == signal_no)
  
  return(data.frame("beta_pbs" = pbs_data |> pull(beta), 
                    "beta_fnf" = fnf_data |> pull(beta),
                    "nompval_pbs" = pbs_data |> pull(nom_pval),
                    "nompval_fnf" = fnf_data |> pull(nom_pval)))
  
}

# Function to determine the order of genotypes in a boxplot from major hom,
# het, to minor hom
determine_geno_order <- function(varID, data){
  
  
  # Get minor allele
  var_minor_allele <- data |> 
    filter(variantID == varID) |> 
    pull(MA) |> 
    unique()
  
  # Get alleles
  alleles <- c(unlist(str_split(varID, ":"))[3], unlist(str_split(varID, ":"))[4])
  
  var_major_allele <- alleles[which(alleles != var_minor_allele)]
  
  # Possible genotypes
  genotypes <- data |> 
    filter(variantID == varID) |> 
    pull(gt_GT_alleles) |> 
    unique()
  
  # Determine het allele order
  possible_hets <- c(paste0(var_major_allele, "/", var_minor_allele),
                     paste0(var_minor_allele, "/", var_major_allele))
  
  var_het <- possible_hets[which(possible_hets %in% genotypes)]
  
  # Make sure order goes from major to minor
  
  
  geno_order <- tibble(variantID = varID,
                       gt_GT_alleles = c(paste0(var_major_allele, "/", var_major_allele),
                                         var_het,
                                         paste0(var_minor_allele, "/", var_minor_allele)),
                       gt_order = c("1", "2", "3"))
  
  return(geno_order)
  
}

# Helper function to convert a (human) gene symbol to its ENSEMBL ID
symbol_to_ensembl <- function(gene_symbol){
  ensembl_id <- AnnotationDbi::select(org.Hs.eg.db, 
                        keys = gene_symbol, 
                        columns = "ENSEMBL", keytype = "SYMBOL") |> pull(ENSEMBL)
}