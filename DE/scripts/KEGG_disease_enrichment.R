library(tidyverse)
library(KEGGREST)
library(org.Hs.eg.db)
library(AnnotationDbi)

# Functions ---------------------------------------------------------------

# Function to grab genes from KEGG disease database entry and convert from
# NCBI to ENSEMBL ID
assignGenestoDisease <- function(diseaseCode){
  dbEntry <- tryCatch({keggGet(diseaseCode)},
                      error = function(cond){
                        return(NULL)
                      })
  if (is.null(dbEntry[[1]]$GENE)){
    return(NA)
  } else {
    genes <- dbEntry[[1]]$GENE
    # Extract hsa number(s)
    gene_ids <- str_extract(genes, "\\[(.*?)\\]")
    # Extract all numbers from string, removing additional letters and symbols
    ncbi_genes <- unlist(str_extract_all(gene_ids, "\\d+"))
    # Convert to ENSEMBL
    ensembl_genes <- tryCatch({AnnotationDbi::select(org.Hs.eg.db, 
                                           keys = ncbi_genes,
                                           keytype = "ENTREZID",
                                           columns = "ENSEMBL")$ENSEMBL},
                              error = function(cond){
                                return(NA)
                              })
    return(ensembl_genes)
    
  }
}

# Function to grab genes from KEGG disease based on pathways in entry
assignGenestoDisease_pathways <- function(diseaseCode){
  
  get_pwGenes <- function(pwid){
    pwEntry <- keggGet(pwid)
    pwGenes <- pwEntry[[1]]$GENE[c(TRUE, FALSE)]
    
    ensembl_genes <- tryCatch({AnnotationDbi::select(org.Hs.eg.db,
                                                     keys = pwGenes,
                                                     keytype = "ENTREZID",
                                                     columns = "ENSEMBL")$ENSEMBL},
                              error = function(cond){
                                return(NA)
                              })
    return(ensembl_genes)
    
  }
  
  

  dbEntry <- tryCatch({keggGet(diseaseCode)},
                      error = function(cond){
                        return(NULL)
                      })
  
  if (is.null(dbEntry[[1]]$PATHWAY)){
    return(NA)
  } else {
    pathways <- gsub("\\([^()]*\\)", "", names(dbEntry[[1]]$PATHWAY))
    print(pathways)
    #pwEntry <- keggGet(pathways)
    #pwGenes <- unlist(pwEntry %>% map("GENE"))[c(TRUE, FALSE)]
    pwGenes <- unlist(lapply(pathways, get_pwGenes))
    
    return(pwGenes)
  }
  
  
}

# Function to perform Wilcoxon rank-sum test for enrichment
performTest <- function(diseaseCode, genes_by_disease, deGenes){
  # Grab genes associated with disease
  diseaseGenes <- genes_by_disease[[diseaseCode]]
  
  # Separate DE genes by which are asssociated with disease and which are not
  genes_in_disease <- intersect(names(deGenes), diseaseGenes)
  scores_in_disease <- deGenes[genes_in_disease]
  
  genes_not_disease <- setdiff(names(deGenes), genes_in_disease)
  scores_not_disease <- deGenes[genes_not_disease]
  
  # Perform Wilcoxon rank-sum test
  if (length(scores_in_disease) > 0){
    wpval <- wilcox.test(scores_in_disease, 
                         scores_not_disease, 
                         alternative = "less")$p.value
  } else {
    wpval <- NA
  }
  
  return(c("wpval" = wpval, 
           "Genes in Term" = length(diseaseGenes),
           "Target Genes in Term" = length(genes_in_disease)))
  
}

# Differential genes ------------------------------------------------------

upsig_de_genes <- read_csv("data/upsig_deGenes_pval01_l2fc2.csv")
upsig_geneList <- upsig_de_genes$padj
names(upsig_geneList) <- upsig_de_genes$gene_id

downsig_de_genes <- read_csv("data/downsig_deGenes_pval01_l2fc2.csv")
downsig_geneList <- downsig_de_genes$padj
names(downsig_geneList) <- downsig_de_genes$gene_id


# Parsing KEGG disease and genes ------------------------------------------

diseaseList <- keggList("disease")
diseaseCodes <- names(diseaseList)

genes_by_disease <- sapply(diseaseCodes, assignGenestoDisease)

# Parsing KEGG disease for pathways and those genes -----------------------

genes_by_pathway_disease <- sapply(diseaseCodes, assignGenestoDisease_pathways)

# Test 1 --------------------------------------------------------------------

upsig_disease <- t(sapply(diseaseCodes, performTest, 
       genes_by_disease = genes_by_disease, deGenes = upsig_geneList)) %>% 
  as.data.frame() %>%
  rownames_to_column(var = "diseaseCode") %>%
  left_join(enframe(diseaseList, name = "diseaseCode", value = "name"), 
            by = "diseaseCode")

upsig_disease_enriched <- upsig_disease %>% filter(wpval < 0.01)

downsig_disease <- t(sapply(diseaseCodes, performTest, 
                          genes_by_disease = genes_by_disease, 
                          deGenes = downsig_geneList)) %>% 
  as.data.frame() %>%
  rownames_to_column(var = "diseaseCode") %>%
  left_join(enframe(diseaseList, name = "diseaseCode", value = "name"), 
            by = "diseaseCode")

downsig_disease_enriched <- downsig_disease %>% filter(wpval < 0.01)


# Test 2 ------------------------------------------------------------------

upsig_disease_pathway <- t(sapply(diseaseCodes, performTest, 
                          genes_by_disease = genes_by_pathway_disease, 
                          deGenes = upsig_geneList)) %>% 
  as.data.frame() %>%
  rownames_to_column(var = "diseaseCode") %>%
  left_join(enframe(diseaseList, name = "diseaseCode", value = "name"), 
            by = "diseaseCode")

upsig_disease_pathway_enriched <- upsig_disease_pathway %>% 
  filter(wpval < 0.01) %>%
  arrange(wpval)

downsig_disease_pathway <- t(sapply(diseaseCodes, performTest, 
                                  genes_by_disease = genes_by_pathway_disease, 
                                  deGenes = downsig_geneList)) %>% 
  as.data.frame() %>%
  rownames_to_column(var = "diseaseCode") %>%
  left_join(enframe(diseaseList, name = "diseaseCode", value = "name"), 
            by = "diseaseCode")

downsig_disease_pathway_enriched <- downsig_disease_pathway %>% 
  filter(wpval < 0.01) %>%
  arrange(wpval)
