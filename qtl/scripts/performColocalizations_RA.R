suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(coloc))

# Input parsing -----------------------------------------------------------
option_list <- list( 
  make_option("--eSNP_eGenes",
              help = "File of significant lead eSNP-eGene pairs with LD information. variantIDs and ld_variantIDs should also have columns named variantID_alpha and ld_variantID_alpha with variantID and ld_variantID alleles put in alphabetical order."),
  make_option(c("--GWAS_window (bp)"), default = 250000, 
              help = "The window size in bp to include on either side of the GWAS lead variant when defining the GWAS signal."),
  make_option("--nominal_prefix",
              help = "The file path and prefix of the nominal QTL results to use for colocalization. These results should be split by chromosome and contain minor allele frequencies for each variant in a column named 'maf'."),
  make_option(c("-p", "--permData"), 
              help = "The file to use to reference QTL permutation pass results for obtaining eGene qvalues."),
  make_option(c("-n", "--sample_size"), 
              help = "The sample size of the QTL study."),
  make_option(c("-o","--outfile"),
              help = "The outfile (.rda) to write results to.")
)

opt <- parse_args(OptionParser(option_list = option_list))

# Functions ---------------------------------------------------------------

parse_gwas_variantID <- function(varRow, type = "lead"){
  
  chrom <- varRow[["chrom"]]
  if (type == "lead"){
    pos <- varRow[["hg38pos"]]
    a1 <- varRow[["EA"]]
    a2 <- varRow[["NEA"]]
  } else {
    pos <- varRow[["ldbuddy_hg38pos"]]
    a1 <- varRow[["ldbuddy_EA"]]
    a2 <- varRow[["ldbuddy_NEA"]]
  }
  
  variantID <- paste0(chrom, ":", pos, ":", paste(sort(c(a1,a2)), collapse = ":"))
  return(variantID)
}

performColocalization <- function(eSNP_eGene, GWAS_leads, QTL_ld, RAtype, 
                                  RApanel, nominalPrefix, qtl_sample_size,
                                  GWAS_window, permData){
  
  # Get signal's eGene id, lead variantID, and signal number
  qtl_variantID_alpha <- eSNP_eGene[["variantID_alpha"]]
  geneid <- eSNP_eGene[["gene_id"]]
  signal_no <- eSNP_eGene[["signal"]]
  
  print(paste0("Performing colocalization for eGene ", geneid, " for signal ", signal_no, "."))
  # Look up qtl_variantID in GWAS_leads ld buddies
  # LD buddies also includes the leads themselves with R2 = 1
  gwas_variant <- GWAS_leads |>  
    filter(ldbuddy_variantID_coloc == qtl_variantID_alpha)
  # Get LD R2 between the GWAS lead variant and the QTL lead variant based on the GWAS LD panel
  gwas_R2 <- gwas_variant$ldbuddy_R2
  
  # Look up lead gwas_variant variantIDs in the QTL lead ld buddies, trying to find
  # both versions of variantIDs
  qtl_R2 <- QTL_ld |>
    # Filter for lead variant/eGene
    filter(variantID_alpha == qtl_variantID_alpha) |> 
    # Find GWAS variant
    filter(ld_variantID_alpha == gwas_variant$variantID_coloc) |> 
    pull(R2)
  
  # Make sure LD is high enough (trying either GWAS or QTL panel)
  if (any(gwas_R2 > 0.5, qtl_R2 > 0.5) & all(!is.na(gwas_R2), !is.na(qtl_R2))){
    
    chr <- gwas_variant$chrom
    
    # GWAS  --------------------------------------------------------------------
    
    # Define GWAS signal region 
    
    min_region <- gwas_variant$hg38pos - GWAS_window
    max_region <- gwas_variant$hg38pos + GWAS_window
      
    # Grab signal from ld buddies
    gwas_sum_chrom <- GWAS_leads |> 
      # Filter for variants within region of lead gwas variant and with
      # valid beta 
      filter(chrom == chr & 
               ldbuddy_hg38pos >= min_region & ldbuddy_hg38pos <= max_region &
               !is.na(ldbuddy_beta)) |>
      # GWAS gives EAF, so convert to MAF
      mutate(MAF = ifelse(ldbuddy_EAF < 0.5, 
                          as.numeric(ldbuddy_EAF), 1 - as.numeric(ldbuddy_EAF))) |> 
      # Determine minor allele based on MAF and EAF
      mutate(MA = ifelse(ldbuddy_EAF == MAF, ldbuddy_EA, ldbuddy_NEA)) |> 
      # Flip sign of beta if NEA is not MA for consistency with QTL data beta (relative to minor allele)
      mutate(BETA_adjusted = ifelse(ldbuddy_NEA != MA, -1*ldbuddy_beta, ldbuddy_beta))
      
    # QTL data -----------------------------------------------------------------
    
    # Pull input nominal QTL data for that independent eGene signal
    qtlData <- read_csv(paste0(nominalPrefix, "_", chr, ".csv")) |> 
      filter(gene_id == geneid & signal == signal_no) |> 
      # Update variantID with alleles in alphabetical order to match GWAS variantIDs
      rowwise() |> 
      mutate(variantID_alpha = paste0(variant_chr, ":", variant_start, ":", 
                                      paste(sort(c(A1, A2)), collapse = ":"))) |> 
      # Square beta standard error to get var beta
      mutate(var_beta = beta_se^2) 
    
    # Sample sizes -----------------------------------------------------------
    
    # Get Case/Control number for gwas
    case_control_sizes <- read_csv("/proj/phanstiel_lab/External/gwas/RA/RA_cases_controls.csv") |> 
      filter(RA == RAtype & ancestry == RApanel)
    fraction_cases <- case_control_sizes$cases/(case_control_sizes$cases + 
                                                      case_control_sizes$controls)
    
    gwasN <- case_control_sizes$cases + case_control_sizes$controls
    
    # QTL
    qtlN <- qtl_sample_size
    
    # run coloc --------------------------------------------------------------
    coloc_result <- tryCatch({coloc.abf(dataset1 = list(pvalues = qtlData$nom_pval,
                                                        N = qtlN,
                                                        MAF = qtlData$maf,
                                                        type = "quant",
                                                        beta = qtlData$beta,
                                                        varbeta = qtlData$var_beta,
                                                        snp = qtlData$variantID_alpha),
                                        dataset2 = list(beta = gwas_sum_chrom$BETA_adjusted,
                                                        s = fraction_cases,
                                                        N = gwasN,
                                                        type = "cc",
                                                        MAF = gwas_sum_chrom$MAF,
                                                        pvalues = gwas_sum_chrom$ldbuddy_pval,
                                                        snp = gwas_sum_chrom$ldbuddy_variantID_coloc))},
                             error = function(cond) {NULL})
    if (!is.null(coloc_result)){
      coloc_result$eGene_id <- geneid
      coloc_result$eGene_name <- eSNP_eGene[["gene_symbol"]]
      coloc_result$signal <- signal_no
      coloc_result$GWAS_lead <- gwas_variant[["rsID"]]
      coloc_result$qtl_variantID <- eSNP_eGene[["variantID"]]
      coloc_result$gwas_R2 <- gwas_R2
      coloc_result$qtl_R2 <- qtl_R2
      
      # Get eGene qvalue from permutation results
      eGene_qval <- read_csv(permData, col_select = c("gene_id", "qval")) |> 
        filter(gene_id == geneid) |> 
        pull(qval)
      
      coloc_result$eGene_qval <- eGene_qval
      
    }    
  } else {
    coloc_result <- NULL
  }
  
  return(coloc_result)
}

# Get independent significant lead eGene-snp pairs to colocalize ---------------

# Read in eSNP_eGenes with LD information 
QTL_ld <- fread(opt$eSNP_eGenes, data.table = FALSE)

print("Read in QTL data")

# Subset for unique lead eSNP_eGene pairs for each independent signal
eGene_snps <- QTL_ld |> 
  dplyr::select(gene_id, gene_symbol, rsID, variantID, variantID_alpha, signal) |> 
  distinct()

# Iterate through GWAS and perform colocalization -------------------------

RAsubtypes <- c("AllRA", "SeropositiveRA")
ra_colocs <- list()

for (subtype in RAsubtypes){
  
  for (panel in c("EAS", "EUR", "multi_ancestry")){
    if (subtype == "AllRA" & panel == "EAS"){
      next
    }
    
    print(paste0("Processing ", panel, " ", subtype, "."))
    
    # Read in leads and LD buddies
    subtype_leads <- read_csv(paste0("/proj/phanstiel_lab/External/gwas/RA/",
                                     subtype, "/", panel, "/leads/",
                                     subtype, "_", panel, "_leads_ld_hg38_final.csv.gz"))
    
    subtype_leads$variantID_coloc <- apply(subtype_leads, 1, 
                                     parse_gwas_variantID, type = "lead")
    subtype_leads$ldbuddy_variantID_coloc <- apply(subtype_leads, 1, 
                                             parse_gwas_variantID, type = "ld")
    
    # Subset lead QTLs for ones found in GWAS subtype leads and LD buddies
    # These are the ones we can test
    eGene_snps_gwasoverlaps <- eGene_snps |> 
      filter(variantID_alpha %in% subtype_leads$ldbuddy_variantID_coloc)
    
    if (nrow(eGene_snps_gwasoverlaps) > 0){
      colocalizations <- apply(eGene_snps_gwasoverlaps, 1, 
                               performColocalization,
                               GWAS_leads = subtype_leads,
                               QTL_ld = QTL_ld,
                               RAtype = subtype,
                               RApanel = panel,
                               nominalPrefix = opt$nominal_prefix,
                               qtl_sample_size = as.numeric(opt$sample_size),
                               GWAS_window = as.numeric(opt$GWAS_window),
                               permData = opt$permData)
      if (!is.null(colocalizations)){
        names(colocalizations) <- eGene_snps_gwasoverlaps$rsID
      }
      
    } else {
      colocalizations <- NULL
    }
    
    ra_colocs[[paste0(subtype, "_", panel)]] <- colocalizations
  }
  

}

save(ra_colocs, file = opt$outfile)

