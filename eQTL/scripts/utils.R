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

#' A function to appropriately read in and name the columns of output
#' from QTLtools cis permutation pass
readQTLtools_perm <- function(filePath){
  
  qtlData <- read_delim(filePath, col_names = FALSE)
  colnames(qtlData) <- c("gene_id", "gene_chr", "gene_start", "gene_end", 
                         "gene_strand", "num_cis_variants", "dist", "variantID",
                         "variant_chr", "variant_start", "variant_end",
                         "dof1", "dof2", "bml1", "bml2", "nom_pval", 
                         "r_squared", "beta", "adj_emp_pval", "adj_beta_pval")
  
  return(qtlData)
}

readQTLtools_nom <- function(filePath){
  qtlData <- read_delim(filePath, col_names = FALSE)
  colnames(qtlData) <- c("gene_id", "gene_chr", "gene_start", "gene_end",
                         "gene_strand", "num_cis_variants", "dist",
                         "variantID", "variant_chr", "variant_start",
                         "variant_end", "nom_pval", "r_squared", "beta",
                         "best_hit")
  return(qtlData)
}

#' A function to convert ENSEMBL gene_ids to SYMBOL with a Bioconductor orgDb 
getGeneSymbol <- function(gene_id, org_db){
  geneSymbols <- AnnotationDbi::select(org_db, keys = gene_id, columns = "SYMBOL", keytype = "ENSEMBL")
  return(geneSymbols)
}

#' A function to calculate the percent overlap of dataset a in dataset b
percentOverlap <- function(a, b){
  intersection <- length(intersect(a,b))
  percentOverlap <- intersection/length(a)
  return(percentOverlap)
}

#' A function to calculate the Jaccard index between dataset a and dataset b
jaccard <- function(a, b) {
  intersection <- length(intersect(a, b))
  union <- length(a) + length(b) - intersection
  jaccard <- intersection/union
  
  return(jaccard)
}

#' Loads an RData file, and returns it
loadRData <- function(fileName){
  
  load(fileName)
  get(ls()[ls() != "fileName"])
}