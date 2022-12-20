library(tidyverse)


# CTL and FNF data covariate options --------------------------------------

conditions <- c("CTL", "FNF")
covariateOptions <- c(paste0(conditions, "_none.txt"),
                      paste0(conditions, "_genoPC_batch.txt"),
                      paste0(conditions, "_genoPC.txt"),
                      paste0(conditions, "_batch.txt"))

allOptions <- read_csv(paste0("output/summary/", covariateOptions)) %>%
  # filter options where discovery was 0
  filter(qval > 0) %>%
  # Add separate columns each describing toggle of covariate options
  mutate(Condition = factor(unlist(lapply(str_split(.$Covariate_Options, "_", n = 2), `[[`, 1)))) %>%
  mutate(RNAmethod = factor(unlist(lapply(str_split(.$Covariate_Options, "_"), `[[`, 2)))) %>%
  mutate(genoPCs = str_detect(.$Covariate_Options, "genoPC")) %>%
  mutate(batch = str_detect(.$Covariate_Options, "batch")) %>%
  mutate(Nk = str_extract(.$Covariate_Options, "k[0-9]")) %>%
  mutate(Name = str_remove(unlist(lapply(str_split(.$Covariate_Options, "_", n = 2), `[[`, 2)), "_FDR")) %>%
  group_by(Condition) %>%
  group_split()


Nks1 <- unique(allOptions[[1]]$Nk)
Nks2 <- unique(allOptions[[2]]$Nk)
common_Nk <- intersect(Nks1, Nks2)


CTL_options <- allOptions[[1]] %>% filter(Nk %in% common_Nk)
FNF_options <- allOptions[[2]] %>% filter(Nk %in% common_Nk)

# GTEx data ---------------------------------------------------------------

# Get list of GTEx tissues
gtex_tissues <- list.files("/pine/scr/n/e/nekramer/GTEx/eqtls/", 
                           pattern = "*.egenes.txt.gz")
gtex_tissues <- gsub(".v8.EUR.egenes.txt.gz", "", gtex_tissues)
  

# Iterate through each GTEx tissue ----------------------------------------
# 
# 
# tissue_percentOverlaps <- list()
# 
# for (tissue in gtex_tissues){
#   # Read in eGenes
#   tissue_eGenes <- read_delim(paste0("/pine/scr/n/e/nekramer/GTEx/eqtls/",
#                                      tissue, ".v8.EUR.egenes.txt.gz"),
#                               delim = "\t") %>%
#     pull(phenotype_id)
#   
#   # Remove decimals from ENS IDs
#   tissue_eGenes <- gsub("\\..*","", tissue_eGenes)
#   
#   
#   # CTL ---------------------------------------------------------------------
# 
#   CTL_percentOverlap <- list()
#   
#   for (option in CTL_options$Name){
#     # Read in eGenes for set of covariate options
#     eGenes_FDR <- read_csv(paste0("output/qtl/CTL_", option, "_perm1Mb_FDR.txt")) %>%
#       filter(FDR <= 0.05) %>%
#       pull(gene_id)
#     
#     eGenes_qval <- read_csv(paste0("output/qtl/CTL_", option, "_perm1Mb_FDR.txt")) %>%
#       filter(qval <= 0.05) %>%
#       pull(gene_id)
#     
#     # Calculate percent overlap
#     percOverlap_FDR <- percentOverlap(eGenes_FDR, tissue_eGenes)
#     percOverlap_qval <- percentOverlap(eGenes_qval, tissue_eGenes)
#     CTL_percentOverlap[[paste0(option, "_FDR")]] <- percOverlap_FDR
#     CTL_percentOverlap[[paste0(option, "_qval")]] <- percOverlap_qval
# 
#   }
#   
# 
#   CTL_percentOverlap <- data.frame("Name" = names(CTL_percentOverlap), 
#                                             "val" = unlist(CTL_percentOverlap))
#   
#   CTL_percentOverlap_qval <- CTL_percentOverlap %>% filter(str_detect(Name, "qval")) %>%
#     dplyr::rename(qval = val) %>%
#     mutate(Name = gsub("_qval", "", .$Name)) %>%
#     mutate(Condition = "CTL")
#   
#   CTL_percentOverlap_BH <- CTL_percentOverlap %>% filter(str_detect(Name, "FDR")) %>%
#     dplyr::rename(FDR = val) %>%
#     mutate(Name = gsub("_FDR", "", .$Name)) %>%
#     mutate(Condition = "CTL")
#   
#   CTL_percentOverlap_all <- left_join(CTL_percentOverlap_qval, 
#                                       CTL_percentOverlap_BH)
#   
#   tissue_percentOverlaps[[paste0(tissue, "_CTL")]] <- CTL_percentOverlap_all
# 
#   
#   
#   # FNF ---------------------------------------------------------------------
#   
#   FNF_percentOverlap <- list()
#   
#   for (option in FNF_options$Name){
#     # Read in eGenes for set of covariate options
#     eGenes_FDR <- read_csv(paste0("output/qtl/FNF_", option, "_perm1Mb_FDR.txt")) %>%
#       filter(FDR <= 0.05) %>%
#       pull(gene_id)
#     
#     eGenes_qval <- read_csv(paste0("output/qtl/FNF_", option, "_perm1Mb_FDR.txt")) %>%
#       filter(qval <= 0.05) %>%
#       pull(gene_id)
#     
#     # Calculate percent overlap
#     percOverlap_FDR <- percentOverlap(eGenes_FDR, tissue_eGenes)
#     percOverlap_qval <- percentOverlap(eGenes_qval, tissue_eGenes)
#     FNF_percentOverlap[[paste0(option, "_FDR")]] <- percOverlap_FDR
#     FNF_percentOverlap[[paste0(option, "_qval")]] <- percOverlap_qval
#     
#   }
#   
#   FNF_percentOverlap <- data.frame("Name" = names(FNF_percentOverlap), 
#                                    "val" = unlist(FNF_percentOverlap))
#   
#   FNF_percentOverlap_qval <- FNF_percentOverlap %>% filter(str_detect(Name, "qval")) %>%
#     dplyr::rename(qval = val) %>%
#     mutate(Name = gsub("_qval", "", .$Name)) %>%
#     mutate(Condition = "FNF")
#   
#   FNF_percentOverlap_BH <- FNF_percentOverlap %>% filter(str_detect(Name, "FDR")) %>%
#     dplyr::rename(FDR = val) %>%
#     mutate(Name = gsub("_FDR", "", .$Name)) %>%
#     mutate(Condition = "FNF")
#   
#   FNF_percentOverlap_all <- left_join(FNF_percentOverlap_qval, 
#                                       FNF_percentOverlap_BH)
#   tissue_percentOverlaps[[paste0(tissue, "_FNF")]] <- FNF_percentOverlap_all
#   
# }
# 
# 
# # Combine and plot --------------------------------------------------------
# 
# percentOverlaps_all <- bind_rows(tissue_percentOverlaps, .id = "Tissue") %>%
#   mutate(Tissue = gsub("[_](CTL|FNF)", "", .$Tissue)) %>%
#   mutate(Tissue = factor(Tissue))
# 
# 
# for (tissue in gtex_tissues){
#   
#   Names_percentOverlaps <- percentOverlaps_all %>% filter(Tissue == tissue) %>% 
#     filter(str_detect(Condition, "CTL")) %>%
#     arrange(FDR) %>%
#     pull(Name)
#   
#   percentOverlaps_data <- percentOverlaps_all %>% filter(Tissue == tissue)
#   percentOverlaps_data$Name <- factor(percentOverlaps_data$Name, levels = Names_percentOverlaps)
#   
#   ggplot(data = percentOverlaps_data) +
#     geom_segment(aes(x = FDR, xend = qval, y = Name, yend = Name), lwd = 0.5) +
#     geom_point(aes(x = FDR, y = Name,color = brewer.pal(n = 3, "Paired")[1]), size = 2) +
#     geom_point(aes(x = qval, y = Name,color = brewer.pal(n = 3, "Paired")[2]), size = 2) +
#     facet_wrap(~Condition,
#                scales = "free_y")+
#     scale_color_identity() +
#     xlab("Percent Overlap") +
#     labs(color = "Legend") +
#     ggtitle(tissue) +
#     theme_minimal() + 
#     theme(axis.title.y=element_blank(),
#           panel.grid.major.y = element_blank(),
#           legend.title = element_blank()) +
#     scale_color_identity(guide = "legend",
#                          labels = c("q-value", "BH"))
#   
#   ggsave(filename = paste0("output/plots/GTEx_", tissue, "_percentOverlap.pdf"))
# }


# ranked p-value splits ---------------------------------------------------

tissue_sliced_percentOverlaps <- list()

for (tissue in gtex_tissues){
  
  # Read in eGenes
  tissue_eGenes <- read_delim(paste0("/pine/scr/n/e/nekramer/GTEx/eqtls/",
                                     tissue, ".v8.EUR.egenes.txt.gz"),
                              delim = "\t") %>% filter(qval <= 0.05) %>%
    pull(phenotype_id)
  
  # Remove decimals from ENS IDs
  tissue_eGenes <- gsub("\\..*","", tissue_eGenes)

  CTL_sliced_percentOverlap <- list()
  for (option in CTL_options$Name){
    # Read in eGenes for set of covariate options and FDR
    eGenes_FDR <- read_csv(paste0("output/qtl/CTL_", option, "_perm1Mb_FDR.txt")) %>%
      filter(FDR <= 0.05)
    n_slices <- ceiling(nrow(eGenes_FDR)/50) 
  
    for (n in seq(1, n_slices, by = 1)){
      n_genes <- n*50
      slice <- eGenes_FDR %>% slice_min(order_by = FDR, n = n_genes) %>%
        pull(gene_id)
    
      # Calculate overlap with tissue
      percOverlap <- percentOverlap(slice, tissue_eGenes)

      CTL_sliced_percentOverlap[[paste0(option, "_top", n_genes, "_FDR")]] <- percOverlap
    }
  
    eGenes_qval <- read_csv(paste0("output/qtl/CTL_", option, "_perm1Mb_FDR.txt")) %>%
      filter(qval <= 0.05)
    n_slices <- ceiling(nrow(eGenes_qval)/50)
  
    for (n in seq(1, n_slices, by = 1)){
      n_genes <- n*50
      slice <- eGenes_qval %>% slice_min(order_by = qval, n = n_genes) %>%
        pull(gene_id)
    
    # Calculate overlap with tissue
    percOverlap <- percentOverlap(slice, tissue_eGenes)

    CTL_sliced_percentOverlap[[paste0(option, "_top", n_genes, "_qval")]] <- percOverlap

    }
  
  }

  CTL_sliced_percentOverlap <- data.frame("Name" = names(CTL_sliced_percentOverlap), 
                                                 "val" = unlist(CTL_sliced_percentOverlap))

  CTL_sliced_percentOverlap_qval <- CTL_sliced_percentOverlap %>% 
    filter(str_detect(Name, "qval")) %>%
    dplyr::rename(qval = val) %>%
    mutate(Name = gsub("_qval", "", .$Name)) %>%
    mutate(NumGenes = str_extract(.$Name, "top[0-9]+")) %>%
    mutate(NumGenes = as.numeric(gsub("top", "", .$NumGenes))) %>%
    mutate(Condition = "CTL")

  CTL_sliced_percentOverlap_BH <- CTL_sliced_percentOverlap %>% 
    filter(str_detect(Name, "FDR")) %>%
    dplyr::rename(FDR = val) %>%
    mutate(Name = gsub("_FDR", "", .$Name)) %>%
    mutate(NumGenes = str_extract(.$Name, "top[0-9]+")) %>%
    mutate(NumGenes = as.numeric(gsub("top", "", .$NumGenes))) %>%
    mutate(Condition = "CTL")

  CTL_sliced_percentOverlap_all <- left_join(CTL_sliced_percentOverlap_qval, CTL_sliced_percentOverlap_BH)
  tissue_sliced_percentOverlaps[[paste0(tissue, "_CTL")]] <- CTL_sliced_percentOverlap_all

  
  FNF_sliced_percentOverlap <- list()
  for (option in FNF_options$Name){
    # Read in eGenes for set of covariate options and FDR
    eGenes_FDR <- read_csv(paste0("output/qtl/FNF_", option, "_perm1Mb_FDR.txt")) %>%
      filter(FDR <= 0.05)
    n_slices <- ceiling(nrow(eGenes_FDR)/50) 
    
    for (n in seq(1, n_slices, by = 1)){
      n_genes <- n*50
      slice <- eGenes_FDR %>% slice_min(order_by = FDR, n = n_genes) %>%
        pull(gene_id)
      
      # Calculate overlap with tissue
      percOverlap <- percentOverlap(slice, tissue_eGenes)
      
      FNF_sliced_percentOverlap[[paste0(option, "_top", n_genes, "_FDR")]] <- percOverlap
    }
    
    eGenes_qval <- read_csv(paste0("output/qtl/FNF_", option, "_perm1Mb_FDR.txt")) %>%
      filter(qval <= 0.05)
    n_slices <- ceiling(nrow(eGenes_qval)/50)
    
    for (n in seq(1, n_slices, by = 1)){
      n_genes <- n*50
      slice <- eGenes_qval %>% slice_min(order_by = qval, n = n_genes) %>%
        pull(gene_id)
      
      # Calculate overlap with tissue
      percOverlap <- percentOverlap(slice, tissue_eGenes)
      
      FNF_sliced_percentOverlap[[paste0(option, "_top", n_genes, "_qval")]] <- percOverlap
      
    }
    
  }
  
  FNF_sliced_percentOverlap <- data.frame("Name" = names(FNF_sliced_percentOverlap), 
                                          "val" = unlist(FNF_sliced_percentOverlap))
  
  FNF_sliced_percentOverlap_qval <- FNF_sliced_percentOverlap %>% 
    filter(str_detect(Name, "qval")) %>%
    dplyr::rename(qval = val) %>%
    mutate(Name = gsub("_qval", "", .$Name)) %>%
    mutate(NumGenes = str_extract(.$Name, "top[0-9]+")) %>%
    mutate(NumGenes = as.numeric(gsub("top", "", .$NumGenes))) %>%
    mutate(Condition = "FNF")
  
  FNF_sliced_percentOverlap_BH <- FNF_sliced_percentOverlap %>% 
    filter(str_detect(Name, "FDR")) %>%
    dplyr::rename(FDR = val) %>%
    mutate(Name = gsub("_FDR", "", .$Name)) %>%
    mutate(NumGenes = str_extract(.$Name, "top[0-9]+")) %>%
    mutate(NumGenes = as.numeric(gsub("top", "", .$NumGenes))) %>%
    mutate(Condition = "FNF")
  
  FNF_sliced_percentOverlap_all <- left_join(FNF_sliced_percentOverlap_qval,
                                             FNF_sliced_percentOverlap_BH)
  tissue_sliced_percentOverlaps[[paste0(tissue, "_FNF")]] <- FNF_sliced_percentOverlap_all

}

saveRDS(tissue_sliced_percentOverlaps, "tissue_sliced_percentOverlaps.rds")

tissue_sliced_percentOverlaps <- read_rds("tissue_sliced_percentOverlaps.rds")

sliced_percentOverlaps_all <- bind_rows(tissue_sliced_percentOverlaps, .id = "Tissue") %>%
  mutate(Tissue = gsub("[_](CTL|FNF)", "", .$Tissue)) %>%
  mutate(Tissue = factor(Tissue)) %>%
  mutate(Name = gsub("_top[0-9]+", "", .$Name))

for (tissue in gtex_tissues){

  percentOverlaps_data <- sliced_percentOverlaps_all %>% filter(Tissue == tissue)

  ggplot(data = percentOverlaps_data, aes(x = NumGenes, y = FDR, group = Name, color = Name)) +
    geom_line() +
    facet_wrap(~Condition,
               scales = "free_y") +
    xlab("Top Genes") +
    ylab("Percent Overlap") +
    labs(color = "Legend") +
    ylim(0, 1) +
    ggtitle(tissue) +
    theme_minimal() +
    theme(panel.grid.major.y = element_blank(),
          legend.title = element_blank())

  ggsave(filename = paste0("output/plots/GTEx_", tissue, "_sliced_percentOverlap.pdf"))
}



# Downsampling GTEx -------------------------------------------------------

# num_eGenes <- data.frame(matrix(nrow = length(gtex_tissues), ncol = 2,
#                                 dimnames = list(NULL, c("Tissue", "num_eGenes"))))

# num_eGenes <- data.frame("Tissue" =  gtex_tissues)
# 
# for (tissue in gtex_tissues){
#   
#   # Read in eGenes
#   num_tissue_eGenes <- read_delim(paste0("/pine/scr/n/e/nekramer/GTEx/eqtls/",
#                                      tissue, ".v8.EUR.egenes.txt.gz"),
#                               delim = "\t") %>% filter(qval <= 0.05) %>%
#   num_eGenes[which(num_eGenes$Tissue == tissue), "num_eGenes"] <- num_tissue_eGenes
# }
# 
# 
# ggplot(data = num_eGenes, aes(x = num_eGenes)) +
#   geom_histogram(bins = 30) +
#   theme_minimal() +
#   theme(axis.title.y = element_blank()) +
#   xlab("Number of eGenes") +
#   ggtitle("GTEx Tissue eGenes")
# 
# ggsave(filename = "output/plots/GTEx_eGenes_histogram.pdf")
