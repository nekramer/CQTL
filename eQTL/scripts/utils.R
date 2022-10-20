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
readQTLtools <- function(filePath){
  
  qtlData <- read_delim(filePath, col_names = FALSE)
  colnames(qtlData) <- c("gene_id", "gene_chr", "gene_start", "gene_end", 
                         "gene_strand", "num_cis_variants", "dist", "variantID",
                         "variant_chr", "variant_start", "variant_end",
                         "dof1", "dof2", "bml1", "bml2", "nom_pval", 
                         "r_squared", "beta", "adj_emp_pval", "adj_beta_val")
  
  return(qtlData)
}
