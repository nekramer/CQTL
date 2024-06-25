suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(coloc))


# Input parsing -----------------------------------------------------------

option_list <- list( 
  make_option("--eSNP_eGenes",
              help = "File of significant lead eSNP-eGene pairs with LD information. variantIDs and ld_variantIDs should also have columns named variantID_alpha and ld_variantID_alpha with variantID and ld_variantID alleles put in alphabetical order."),
  make_option(c("-g", "--gwas"), default = "/work/users/n/e/nekramer/External/gwas/Boer_reprocessed/",
              help = "The path to OA GWAS data to use for colocalization. This path should have subfolders for each OA subtype."),
  make_option(c("-l", "--ld_panel"), default = "EUR", 
              help = "A string, either 'EUR' (European) or 'ALL' (All), indicating which 1000G LD panel to use with the GWAS data."),
  make_option(c("--GWAS_signal_thresholding"), default = "region",
              help = "A string, either 'region' (indicating genomic base pair window) or 'ld' (indicating ld thresholding), defining how to subset GWAS variants to define the signal to be colocalized."),
  make_option(c("--GWAS_window (bp)"), default = 1000000, 
              help = "If --GWAS_signal_thresholding='region', the window size in bp to include on either side of the GWAS lead variant when defining the GWAS signal."),
  make_option(c("--GWAS_ld (R2)"), default = 0.6,
              help = "If --GWAS_signal_thresholding='ld', the ld R2 threshold to use for including variants relative to the lead GWAS variant when defining the GWAS signal."),
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
  
  if (type == "lead"){
    chr_pos <- varRow[["CHR:hg38POS"]]
    a1 <- varRow[["EA"]]
    a2 <- varRow[["NEA"]]
  } else {
    chr_pos <- varRow[["ldbuddy_CHR:hg38POS"]]
    a1 <- varRow[["ldbuddy_ref"]]
    a2 <- varRow[["ldbuddy_alt"]]
  }
  
  variantID <- paste0("chr", chr_pos, ":", paste(sort(c(a1,a2)), collapse = ":"))
  return(variantID)
}

performColocalization <- function(eSNP_eGene, GWAS_leads, QTL_ld, OAtype, 
                                  nominalPrefix, gwasPath, qtl_sample_size,
                                  GWAS_signal_thresholding, GWAS_window, 
                                  GWAS_ld, permData){
  
  # Get signal's eGene id, lead variantID, and signal number
  qtl_variantID_alpha <- eSNP_eGene[["variantID_alpha"]]
  geneid <- eSNP_eGene[["gene_id"]]
  signal_no <- eSNP_eGene[["signal"]]
  
  print(paste0("Performing colocalization for eGene ", geneid, " for signal ", signal_no, "."))
  # Look up qtl_variantID in GWAS_leads ld buddies
  # LD buddies also includes the leads themselves with R2 = 1
  gwas_variant <- GWAS_leads |>  
    filter(ldbuddy_variantID == qtl_variantID_alpha)
  # Get LD R2 between the GWAS lead variant and the QTL lead variant based on the GWAS LD panel
  gwas_R2 <- gwas_variant$ldbuddy_R2
  
  # Look up lead gwas_variant variantIDs in the QTL lead ld buddies, trying to find
  # both versions of variantIDs
  qtl_R2 <- QTL_ld |>
    # Filter for lead variant/eGene
    filter(variantID_alpha == qtl_variantID_alpha) |> 
    # Find GWAS variant
    filter(ld_variantID_alpha == gwas_variant$variantID) |> 
    pull(R2)
  
  # Make sure LD is high enough (trying either GWAS or QTL panel)
  if (any(gwas_R2 > 0.5, qtl_R2 > 0.5) & all(!is.na(gwas_R2), !is.na(qtl_R2))){
  
    chrom <- gwas_variant$chrom
    
    # GWAS  --------------------------------------------------------------------
    
    # Define GWAS signal region based on user input
    if (GWAS_signal_thresholding == "region"){
      min_region <- gwas_variant$hg38pos - GWAS_window
      max_region <- gwas_variant$hg38pos + GWAS_window
      
      # Grab gwas sumstats for that chromosome and OA type 
      gwas_sum_chrom <- fread(paste0(gwasPath,
                                     OAtype,
                                     "/summary_stats/",
                                     OAtype, "_chr",
                                     chrom, ".csv"), data.table = FALSE) |> 
        # Filter for variants within region of lead gwas variant
        filter(hg38pos >= min_region & hg38pos <= max_region) |>
        # GWAS gives EAF, so convert to MAF
        mutate(MAF = ifelse(EAF < 0.5, 
                            EAF, 1 - EAF)) |> 
        # Determine minor allele based on MAF and EAF
        mutate(MA = ifelse(EAF == MAF, EA, NEA)) |> 
        # Flip sign of beta if NEA is not MA for consistency with QTL data beta (relative to minor allele)
        mutate(BETA_adjusted = ifelse(NEA != MA, -1*BETA, BETA)) |>
        # Add variantID column with position and alphabetical alleles to match
        # QTL variantIDs
        rowwise() |> 
        mutate(variantID = paste0("chr", `CHR:hg38POS`, ":", 
                                  paste(sort(c(EA, NEA)), collapse = ":")))
      
    } else if (GWAS_signal_thresholding == "ld"){
      
      # Get ld buddies based on GWAS lead
      gwas_variant_LD <- GWAS_leads |> 
        filter(variantID == gwas_variant[["variantID"]]) |> 
        filter(ldbuddy_R2 > GWAS_ld)
      
      # Grab gwas sumstats for that chromosome and OA type 
      gwas_sum_chrom <- fread(paste0(gwasPath,
                                     OAtype,
                                     "/summary_stats/",
                                     OAtype, "_chr",
                                     chrom, ".csv"), data.table = FALSE) |> 
        # filter for variants based on LD
        filter(`CHR:hg38POS` %in% gwas_variant_LD$`ldbuddy_CHR:hg38POS` & !is.na(hg38pos)) |> 
        # GWAS gives EAF, so convert to MAF
        mutate(MAF = ifelse(EAF < 0.5, 
                            EAF, 1 - EAF)) |> 
        # Determine minor allele based on MAF and EAF
        mutate(MA = ifelse(EAF == MAF, EA, NEA)) |> 
        # Flip sign of beta if NEA is not MA for consistency with QTL data beta (relative to minor allele)
        mutate(BETA_adjusted = ifelse(NEA != MA, -1*BETA, BETA)) |>
        # Add variantID column with position and alphabetical alleles to match
        # QTL variantIDs
        rowwise() |> 
        mutate(variantID = paste0("chr", `CHR:hg38POS`, ":", 
                                  paste(sort(c(EA, NEA)), collapse = ":")))
      
    }
    
    # QTL data -----------------------------------------------------------------
    
    # Pull input nominal QTL data for that independent eGene signal
    qtlData <- read_csv(paste0(nominalPrefix, "_chr", chrom, ".csv")) |> 
      filter(gene_id == geneid & signal == signal_no) |> 
      # Update variantID with alleles in alphabetical order to match GWAS variantIDs
      rowwise() |> 
      mutate(variantID_alpha = paste0(variant_chr, ":", variant_start, ":", 
                                      paste(sort(c(A1, A2)), collapse = ":"))) |> 
      # Square beta standard error to get var beta
      mutate(var_beta = beta_se^2) 
    
    # Sample sizes -----------------------------------------------------------
    
    # Get Case/Control number for gwas
    case_control_sizes <- read_csv(paste0(gwasPath, "Case_Control_sampleSizes.csv")) |> 
      filter(OAsubtype == OAtype)
    fraction_cases <- case_control_sizes$Max_Cases/(case_control_sizes$Max_Cases + 
                                                      case_control_sizes$Max_Controls)
    
    gwasN <- case_control_sizes$Max_Cases + case_control_sizes$Max_Controls
    
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
                                               pvalues = gwas_sum_chrom$p,
                                               snp = gwas_sum_chrom$variantID))},
                             error = function(cond) {NULL})
      if (!is.null(coloc_result)){
        coloc_result$eGene_id <- geneid
        coloc_result$eGene_name <- eSNP_eGene[["gene_symbol"]]
        coloc_result$signal <- signal_no
        coloc_result$GWAS_lead <- gwas_variant[["rsID"]]
        coloc_result$GWAS_lead_variantID <- gwas_variant[["variantID"]]
        coloc_result$qtl_variantID <- eSNP_eGene[["variantID"]]
        
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

OAsubtypes <- c("AllOA", "FingerOA", "HandOA", "HipOA", "KneeHipOA", "KneeOA",
                "THR", "ThumbOA", "TJR", "TKR")
all_colocs <- list()

for (subtype in OAsubtypes){
  print(paste0("Processing ", subtype))
  
  # Read in corresponding OA subtype leads and their LD buddies
  # based on input GWAS directory and LD panel
  subtype_leads <- read_csv(paste0(opt$gwas, subtype, "/leads/", opt$ld_panel, "_", subtype, 
                                   "_leads_ld_final.csv"),
    col_types = "cdddccccddddddddddddddddddddddddddddddddddddddddddccdcccd")
  
  subtype_leads$variantID <- apply(subtype_leads, 1, 
                                           parse_gwas_variantID, type = "lead")
  subtype_leads$ldbuddy_variantID <- apply(subtype_leads, 1, 
                                 parse_gwas_variantID, type = "ld")

  
  # Subset lead QTLs for ones found in GWAS subtype leads and LD buddies
  # These are the ones we can test
  eGene_snps_gwasoverlaps <- eGene_snps |> 
    filter(variantID_alpha %in% subtype_leads$ldbuddy_variantID)
  
  
  if (nrow(eGene_snps_gwasoverlaps) > 0){
    colocalizations <- apply(eGene_snps_gwasoverlaps, 1, 
                             performColocalization,
                             GWAS_leads = subtype_leads,
                             QTL_ld = QTL_ld,
                             OAtype = subtype,
                             nominalPrefix = opt$nominal_prefix,
                             gwasPath = opt$gwas,
                             qtl_sample_size = as.numeric(opt$sample_size),
                             GWAS_signal_thresholding = opt$GWAS_signal_thresholding,
                             GWAS_window = as.numeric(opt$GWAS_window),
                             GWAS_ld = as.numeric(opt$GWAS_ld),
                             permData = opt$permData)
    if (!is.null(colocalizations)){
      names(colocalizations) <- eGene_snps_gwasoverlaps$rsID
    }
    
  } else {
    colocalizations <- NULL
  }
  
  all_colocs[[subtype]] <- colocalizations
}

save(all_colocs, file = opt$outfile)
