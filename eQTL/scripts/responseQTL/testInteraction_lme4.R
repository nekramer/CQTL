#' Test for interaction between genotype and condition using ANOVA and lme4
#'
#' @param gene_id Tested gene id
#' @param variant_id Tested SNP id
#' @param expression_matrix expression matrix
#' @param sample_metadata data frame with sample metadata and covariates
#' @param geno_matrix genotype matrix from read.vcfR
#' @param qtl_formula Formula for the model with just genotype and condition terms
#' @param interaction_formula Formula for the model with interaction term between genotype and condition
#'
#' @return Returns a list containing the results, model, and p-value
testInteractionLme4 <- function(gene_id, 
                                variant_id, 
                                expression_matrix, 
                                sample_metadata, 
                                geno_matrix, 
                                qtl_formula,
                                interaction_formula){
  
  # Grab expression of gene_id from expression_matrix
  expression_data <- expression_matrix[gene_id,]
  expression_data <- data.frame(Sample = names(expression_data),
                                expression = expression_data)
  
  # Get genotype information of snp_id from vcf_data
  geno_data <- geno_matrix[variant_id,]
  geno_data <- data.frame(Donor = names(geno_data),
                          genotype = geno_data)
  
  # Join expression_data and genotype_data with sample_metadata
  model_data <- sample_metadata %>%
    left_join(expression_data) %>%
    left_join(geno_data)
  
  # Apply two models to the data and compare them using anova
  no_interaction <- lme4::lmer(qtl_formula, model_data, REML = FALSE)
  interaction <- lme4::lmer(interaction_formula, model_data, REML = FALSE)
  result <- anova(no_interaction, interaction)
  
  # Return results in a list
  return(list(result = result,
              qtl_model = no_interaction,
              interaction = interaction,
              pval = result[[8]][2]))
}