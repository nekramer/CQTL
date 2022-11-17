library(tidyverse)


# FUNCTIONS ---------------------------------------------------------------

# Calculate Jaccard index
jaccard <- function(a, b) {
  intersection <- length(intersect(a, b))
  union <- length(a) + length(b) - intersection
  jaccard <- intersection/union
  
  return(jaccard)
}

# Calculate a normalized Jaccard index?
jaccard_norm <- function(a, b) {
  intersection <- length(intersect(a, b))
  union <- length(a) + length(b) - intersection
  jaccard <- intersection/union
  
  max_jaccard <- min(length(a), length(b))/union
  norm_jaccard <- jaccard/max_jaccard
  
  return(norm_jaccard)
}



# Calculate total percent overlap of our data
percentOverlap <- function(a, b){
  intersection <- length(intersect(a,b))
  percentOverlap <- intersection/length(a)
  return(percentOverlap)
}

# Read in Steinberg et al low grade and high grade eQTLs, getting distinct 
# eGenes
lowgrade_eGenes <- read_delim("/proj/phanstiel_lab/External/public/Steinberg_2021/eQTL_LowGrade_FastQTL.txt", 
                        delim = "\t", 
                        col_types = "ccddddddddcdddddddddcc") %>% 
  distinct(Gene, .keep_all = TRUE) %>%
  separate(Gene, into = c("gene_name", "gene_id"), sep = "_")
highgrade_eGenes <- read_delim("/proj/phanstiel_lab/External/public/Steinberg_2021/eQTL_HighGrade_FastQTL.txt", 
                         delim = "\t", 
                         col_types = "ccddddddddcdddddddddcc") %>%
  distinct(Gene, .keep_all = TRUE) %>%
  separate(Gene, into = c("gene_name", "gene_id"), sep = "_")

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

# CTL calculations --------------------------------------------------------
CTL_jaccards_lowgrade <- list()
CTL_jaccards_highgrade <- list()

CTL_percentOverlap_lowgrade <- list()
CTL_percentOverlap_highgrade <- list()

for (option in CTL_options$Name){
  # Read in eGenes for set of covariate options
  eGenes_FDR <- read_csv(paste0("output/qtl/CTL_", option, "_perm1Mb_FDR.txt")) %>%
    filter(FDR <= 0.05) %>%
    pull(gene_name)
  
  eGenes_qval <- read_csv(paste0("output/qtl/CTL_", option, "_perm1Mb_FDR.txt")) %>%
    filter(qval <= 0.05) %>%
    pull(gene_name)
  
  # Calculate jaccard index with low-grade eQTLs
  jaccard_FDR_lowgrade <- jaccard(eGenes_FDR, lowgrade_eGenes$gene_name)
  jaccard_qval_lowgrade <- jaccard(eGenes_qval, lowgrade_eGenes$gene_name)
  # jaccard_FDR_lowgrade <- jaccard_norm(eGenes_FDR, lowgrade_eGenes$gene_name)
  # jaccard_qval_lowgrade <- jaccard_norm(eGenes_qval, lowgrade_eGenes$gene_name)
  CTL_jaccards_lowgrade[[paste0(option, "_FDR")]] <- jaccard_FDR_lowgrade
  CTL_jaccards_lowgrade[[paste0(option, "_qval")]] <- jaccard_qval_lowgrade
  
  
  # Calculate percent overlap with low-grade eQTLs
  percOverlap_FDR_lowgrade <- percentOverlap(eGenes_FDR, lowgrade_eGenes$gene_name)
  percOverlap_qval_lowgrade <- percentOverlap(eGenes_qval, lowgrade_eGenes$gene_name)
  CTL_percentOverlap_lowgrade[[paste0(option, "_FDR")]] <- percOverlap_FDR_lowgrade
  CTL_percentOverlap_lowgrade[[paste0(option, "_qval")]] <- percOverlap_qval_lowgrade
  
  # Calculate jaccard index with high-grade eQTLs
  jaccard_FDR_highgrade <- jaccard(eGenes_FDR, highgrade_eGenes$gene_name)
  jaccard_qval_highgrade <- jaccard(eGenes_qval, highgrade_eGenes$gene_name)
  # jaccard_FDR_highgrade <- jaccard_norm(eGenes_FDR, highgrade_eGenes$gene_name)
  # jaccard_qval_highgrade <- jaccard_norm(eGenes_qval, highgrade_eGenes$gene_name)
  CTL_jaccards_highgrade[[paste0(option, "_FDR")]] <- jaccard_FDR_highgrade
  CTL_jaccards_highgrade[[paste0(option, "_qval")]] <- jaccard_qval_highgrade
  
  
  # Calculate percent overlap with high-grade eQTLs
  percOverlap_FDR_highgrade <- percentOverlap(eGenes_FDR, highgrade_eGenes$gene_name)
  percOverlap_qval_highgrade <- percentOverlap(eGenes_qval, highgrade_eGenes$gene_name)
  CTL_percentOverlap_highgrade[[paste0(option, "_FDR")]] <- percOverlap_FDR_highgrade
  CTL_percentOverlap_highgrade[[paste0(option, "_qval")]] <- percOverlap_qval_highgrade
}

CTL_jaccards_lowgrade <- data.frame("Name" = names(CTL_jaccards_lowgrade), 
                           "val" = unlist(CTL_jaccards_lowgrade))
  
CTL_jaccards_lowgrade_qval <- CTL_jaccards_lowgrade %>% filter(str_detect(Name, "qval")) %>%
  dplyr::rename(qval = val) %>%
  mutate(Name = gsub("_qval", "", .$Name)) %>%
  mutate(Condition = "CTL") %>%
  mutate(Grade = "low")

CTL_jaccards_lowgrade_BH <- CTL_jaccards_lowgrade %>% filter(str_detect(Name, "FDR")) %>%
  dplyr::rename(FDR = val) %>%
  mutate(Name = gsub("_FDR", "", .$Name)) %>%
  mutate(Condition = "CTL") %>%
  mutate(Grade = "low")

CTL_jaccards_highgrade <- data.frame("Name" = names(CTL_jaccards_highgrade), 
                                    "val" = unlist(CTL_jaccards_highgrade))

CTL_jaccards_highgrade_qval <- CTL_jaccards_highgrade %>% filter(str_detect(Name, "qval")) %>%
  dplyr::rename(qval = val) %>%
  mutate(Name = gsub("_qval", "", .$Name)) %>%
  mutate(Condition = "CTL") %>%
  mutate(Grade = "high")

CTL_jaccards_highgrade_BH <- CTL_jaccards_highgrade %>% filter(str_detect(Name, "FDR")) %>%
  dplyr::rename(FDR = val) %>%
  mutate(Name = gsub("_FDR", "", .$Name)) %>%
  mutate(Condition = "CTL") %>%
  mutate(Grade = "high")


CTL_jaccards_all <- bind_rows(left_join(CTL_jaccards_lowgrade_qval,
                                        CTL_jaccards_lowgrade_BH) %>% arrange(FDR),
                              left_join(CTL_jaccards_highgrade_qval,
                                        CTL_jaccards_highgrade_BH) %>% arrange(FDR))

CTL_jaccards_all$Name <- factor(CTL_jaccards_all$Name, 
                                levels = unique(CTL_jaccards_all$Name))



CTL_percentOverlap_lowgrade <- data.frame("Name" = names(CTL_percentOverlap_lowgrade), 
                                    "val" = unlist(CTL_percentOverlap_lowgrade))

CTL_percentOverlap_lowgrade_qval <- CTL_percentOverlap_lowgrade %>% filter(str_detect(Name, "qval")) %>%
  dplyr::rename(qval = val) %>%
  mutate(Name = gsub("_qval", "", .$Name)) %>%
  mutate(Condition = "CTL") %>%
  mutate(Grade = "low")

CTL_percentOverlap_lowgrade_BH <- CTL_percentOverlap_lowgrade %>% filter(str_detect(Name, "FDR")) %>%
  dplyr::rename(FDR = val) %>%
  mutate(Name = gsub("_FDR", "", .$Name)) %>%
  mutate(Condition = "CTL") %>%
  mutate(Grade = "low")

CTL_percentOverlap_highgrade <- data.frame("Name" = names(CTL_percentOverlap_highgrade), 
                                     "val" = unlist(CTL_percentOverlap_highgrade))

CTL_percentOverlap_highgrade_qval <- CTL_percentOverlap_highgrade %>% filter(str_detect(Name, "qval")) %>%
  dplyr::rename(qval = val) %>%
  mutate(Name = gsub("_qval", "", .$Name)) %>%
  mutate(Condition = "CTL") %>%
  mutate(Grade = "high")

CTL_percentOverlap_highgrade_BH <- CTL_percentOverlap_highgrade %>% filter(str_detect(Name, "FDR")) %>%
  dplyr::rename(FDR = val) %>%
  mutate(Name = gsub("_FDR", "", .$Name)) %>%
  mutate(Condition = "CTL") %>%
  mutate(Grade = "high")


CTL_percentOverlap_all <- bind_rows(left_join(CTL_percentOverlap_lowgrade_qval,
                                        CTL_percentOverlap_lowgrade_BH) %>% arrange(FDR),
                              left_join(CTL_percentOverlap_highgrade_qval,
                                        CTL_percentOverlap_highgrade_BH) %>% arrange(FDR))

CTL_percentOverlap_all$Name <- factor(CTL_percentOverlap_all$Name, 
                                levels = unique(CTL_percentOverlap_all$Name))

# FNF calculations --------------------------------------------------------


FNF_jaccards_lowgrade <- list()
FNF_jaccards_highgrade <- list()

FNF_percentOverlap_lowgrade <- list()
FNF_percentOverlap_highgrade <- list()

for (option in FNF_options$Name){
  # Read in eGenes for set of covariate options
  eGenes_FDR <- read_csv(paste0("output/qtl/FNF_", option, "_perm1Mb_FDR.txt")) %>%
    filter(FDR <= 0.05) %>%
    pull(gene_name)
  
  eGenes_qval <- read_csv(paste0("output/qtl/FNF_", option, "_perm1Mb_FDR.txt")) %>%
    filter(qval <= 0.05) %>%
    pull(gene_name)
  
  # Calculate jaccard index with high-grade eQTLs
  jaccard_FDR_highgrade <- jaccard(eGenes_FDR, highgrade_eGenes$gene_name)
  jaccard_qval_highgrade <- jaccard(eGenes_qval, highgrade_eGenes$gene_name)
  
  # jaccard_FDR_highgrade <- jaccard_norm(eGenes_FDR, highgrade_eGenes$gene_name)
  # jaccard_qval_highgrade <- jaccard_norm(eGenes_qval, highgrade_eGenes$gene_name)
  
  FNF_jaccards_highgrade[[paste0(option, "_FDR")]] <- jaccard_FDR_highgrade
  FNF_jaccards_highgrade[[paste0(option, "_qval")]] <- jaccard_qval_highgrade
  
  # Calculate percent overlap with high-grade eQTLs
  percOverlap_FDR_highgrade <- percentOverlap(eGenes_FDR, highgrade_eGenes$gene_name)
  percOverlap_qval_lowgrade <- percentOverlap(eGenes_qval, highgrade_eGenes$gene_name)
  FNF_percentOverlap_highgrade[[paste0(option, "_FDR")]] <- percOverlap_FDR_highgrade
  FNF_percentOverlap_highgrade[[paste0(option, "_qval")]] <- percOverlap_qval_highgrade
  
  
  # Calculate jaccard index with low-grade eQTLs
  jaccard_FDR_lowgrade <- jaccard(eGenes_FDR, lowgrade_eGenes$gene_name)
  jaccard_qval_lowgrade <- jaccard(eGenes_qval, lowgrade_eGenes$gene_name)
  # jaccard_FDR_lowgrade <- jaccard_norm(eGenes_FDR, lowgrade_eGenes$gene_name)
  # jaccard_qval_lowgrade <- jaccard_norm(eGenes_qval, lowgrade_eGenes$gene_name)
  FNF_jaccards_lowgrade[[paste0(option, "_FDR")]] <- jaccard_FDR_lowgrade
  FNF_jaccards_lowgrade[[paste0(option, "_qval")]] <- jaccard_qval_lowgrade
  
  
  # Calculate percent overlap with low-grade eQTLs
  percOverlap_FDR_lowgrade <- percentOverlap(eGenes_FDR, lowgrade_eGenes$gene_name)
  percOverlap_qval_lowgrade <- percentOverlap(eGenes_qval, lowgrade_eGenes$gene_name)
  FNF_percentOverlap_lowgrade[[paste0(option, "_FDR")]] <- percOverlap_FDR_lowgrade
  FNF_percentOverlap_lowgrade[[paste0(option, "_qval")]] <- percOverlap_qval_lowgrade
}

FNF_jaccards_lowgrade <- data.frame("Name" = names(FNF_jaccards_lowgrade), 
                           "val" = unlist(FNF_jaccards_lowgrade))
FNF_jaccards_highgrade <- data.frame("Name" = names(FNF_jaccards_highgrade), 
                                    "val" = unlist(FNF_jaccards_highgrade))

FNF_jaccards_lowgrade_qval <- FNF_jaccards_lowgrade %>% filter(str_detect(Name, "qval")) %>%
  dplyr::rename(qval = val) %>%
  mutate(Name = gsub("_qval", "", .$Name)) %>%
  mutate(Condition = "FNF") %>%
  mutate(Grade = "low")

FNF_jaccards_lowgrade_BH <- FNF_jaccards_lowgrade %>% filter(str_detect(Name, "FDR")) %>%
  dplyr::rename(FDR = val) %>%
  mutate(Name = gsub("_FDR", "", .$Name)) %>%
  mutate(Condition = "FNF") %>%
  mutate(Grade = "low")

FNF_jaccards_highgrade_qval <- FNF_jaccards_highgrade %>% 
  filter(str_detect(Name, "qval")) %>%
  dplyr::rename(qval = val) %>%
  mutate(Name = gsub("_qval", "", .$Name)) %>%
  mutate(Condition = "FNF") %>%
  mutate(Grade = "high")

FNF_jaccards_highgrade_BH <- FNF_jaccards_highgrade %>% 
  filter(str_detect(Name, "FDR")) %>%
  dplyr::rename(FDR = val) %>%
  mutate(Name = gsub("_FDR", "", .$Name)) %>%
  mutate(Condition = "FNF") %>%
  mutate(Grade = "high")

FNF_jaccards_all <- bind_rows(left_join(FNF_jaccards_lowgrade_qval,
                                        FNF_jaccards_lowgrade_BH),
                              left_join(FNF_jaccards_highgrade_qval,
                                        FNF_jaccards_highgrade_BH))
FNF_jaccards_all$Name <- factor(FNF_jaccards_all$Name, 
                                levels = unique(CTL_jaccards_all$Name))

FNF_percentOverlap_lowgrade <- data.frame("Name" = names(FNF_percentOverlap_lowgrade), 
                                    "val" = unlist(FNF_percentOverlap_lowgrade))
FNF_percentOverlap_highgrade <- data.frame("Name" = names(FNF_percentOverlap_highgrade), 
                                     "val" = unlist(FNF_percentOverlap_highgrade))

FNF_percentOverlap_lowgrade_qval <- FNF_percentOverlap_lowgrade %>% filter(str_detect(Name, "qval")) %>%
  dplyr::rename(qval = val) %>%
  mutate(Name = gsub("_qval", "", .$Name)) %>%
  mutate(Condition = "FNF") %>%
  mutate(Grade = "low")

FNF_percentOverlap_lowgrade_BH <- FNF_percentOverlap_lowgrade %>% filter(str_detect(Name, "FDR")) %>%
  dplyr::rename(FDR = val) %>%
  mutate(Name = gsub("_FDR", "", .$Name)) %>%
  mutate(Condition = "FNF") %>%
  mutate(Grade = "low")

FNF_percentOverlap_highgrade_qval <- FNF_percentOverlap_highgrade %>% 
  filter(str_detect(Name, "qval")) %>%
  dplyr::rename(qval = val) %>%
  mutate(Name = gsub("_qval", "", .$Name)) %>%
  mutate(Condition = "FNF") %>%
  mutate(Grade = "high")

FNF_percentOverlap_highgrade_BH <- FNF_percentOverlap_highgrade %>% 
  filter(str_detect(Name, "FDR")) %>%
  dplyr::rename(FDR = val) %>%
  mutate(Name = gsub("_FDR", "", .$Name)) %>%
  mutate(Condition = "FNF") %>%
  mutate(Grade = "high")

FNF_percentOverlap_all <- bind_rows(left_join(FNF_percentOverlap_lowgrade_qval,
                                        FNF_percentOverlap_lowgrade_BH),
                              left_join(FNF_percentOverlap_highgrade_qval,
                                        FNF_percentOverlap_highgrade_BH))
FNF_percentOverlap_all$Name <- factor(FNF_percentOverlap_all$Name, 
                                levels = unique(CTL_percentOverlap_all$Name))

# Joining and plotting Jaccards -------------------------------------------

jaccards_all <- bind_rows(CTL_jaccards_all, FNF_jaccards_all) %>%
  mutate(Grade = factor(Grade, levels = c("low", "high")))

ggplot(data = jaccards_all) +
  geom_segment(aes(x = FDR, xend = qval, y = Name, yend = Name), lwd = 0.5) +
  geom_point(aes(x = FDR, y = Name,color = brewer.pal(n = 3, "Paired")[1]), size = 2) +
  geom_point(aes(x = qval, y = Name,color = brewer.pal(n = 3, "Paired")[2]), size = 2) +
  facet_wrap(~Condition + Grade,
             scales = "free_y")+
  scale_color_identity() +
  xlab("Jaccard Index") +
  labs(color = "Legend") +
  theme_minimal() + 
  theme(axis.title.y=element_blank(),
        panel.grid.major.y = element_blank(),
        legend.title = element_blank()) +
  scale_color_identity(guide = "legend",
                       labels = c("q-value", "BH"))
ggsave(filename = "output/plots/Steinberg_jaccard_allcombos.pdf")  


# Joining and plotting percent overlaps -----------------------------------
percentOverlap_all <- bind_rows(CTL_percentOverlap_all, FNF_percentOverlap_all) %>%
  mutate(Grade = factor(Grade, levels = c("low", "high")))

ggplot(data = percentOverlap_all) +
  geom_segment(aes(x = FDR, xend = qval, y = Name, yend = Name), lwd = 0.5) +
  geom_point(aes(x = FDR, y = Name, color = brewer.pal(n = 3, "Paired")[1]), size = 2) +
  geom_point(aes(x = qval, y = Name, color = brewer.pal(n = 3, "Paired")[2]), size = 2) +
  facet_wrap(~Condition + Grade,
             scales = "free_y")+
  scale_color_identity() +
  xlab("Percent Overlap") +
  labs(color = "Legend") +
  theme_minimal() + 
  theme(axis.title.y=element_blank(),
        panel.grid.major.y = element_blank(),
        legend.title = element_blank()) +
  scale_color_identity(guide = "legend",
                       labels = c("q-value", "BH"))
ggsave(filename = "output/plots/Steinberg_percOverlap_allcombos.pdf") 


# Ranking by p-value ------------------------------------------------------

CTL_sliced_percentOverlap_lowgrade <- list()
CTL_sliced_percentOverlap_highgrade <- list()
for (option in CTL_options$Name){
  # Read in eGenes for set of covariate options and FDR
  eGenes_FDR <- read_csv(paste0("output/qtl/CTL_", option, "_perm1Mb_FDR.txt")) %>%
    filter(FDR <= 0.05)
  n_slices <- ceiling(nrow(eGenes_FDR)/50) 
  
  for (n in seq(1, n_slices, by = 1)){
    n_genes <- n*50
    slice <- eGenes_FDR %>% slice_min(order_by = FDR, n = n_genes) %>%
      pull(gene_name)
    
    # Calculate overlap with low grade
    percOverlap_lowgrade <- percentOverlap(slice, lowgrade_eGenes$gene_name)
    # Calculate overlap with high grade
    percOverlap_highgrade <- percentOverlap(slice, highgrade_eGenes$gene_name)
    CTL_sliced_percentOverlap_lowgrade[[paste0(option, "_top", n_genes, "_FDR")]] <- percOverlap_lowgrade
    CTL_sliced_percentOverlap_highgrade[[paste0(option, "_top", n_genes, "_FDR")]] <- percOverlap_highgrade
  }
  
  eGenes_qval <- read_csv(paste0("output/qtl/CTL_", option, "_perm1Mb_FDR.txt")) %>%
    filter(qval <= 0.05)
  n_slices <- ceiling(nrow(eGenes_qval)/50)
  
  for (n in seq(1, n_slices, by = 1)){
    n_genes <- n*50
    slice <- eGenes_qval %>% slice_min(order_by = qval, n = n_genes) %>%
      pull(gene_name)
    
    # Calculate overlap with low grade
    percOverlap_lowgrade <- percentOverlap(slice, lowgrade_eGenes$gene_name)
    # Calculate overlap with high grade
    percOverlap_highgrade <- percentOverlap(slice, highgrade_eGenes$gene_name)
    CTL_sliced_percentOverlap_lowgrade[[paste0(option, "_top", n_genes, "_qval")]] <- percOverlap_lowgrade
    CTL_sliced_percentOverlap_highgrade[[paste0(option, "_top", n_genes, "_qval")]] <- percOverlap_highgrade
  }
  
}

CTL_sliced_percentOverlap_lowgrade <- data.frame("Name" = names(CTL_sliced_percentOverlap_lowgrade), 
                                                 "val" = unlist(CTL_sliced_percentOverlap_lowgrade))

CTL_sliced_percentOverlap_lowgrade_qval <- CTL_sliced_percentOverlap_lowgrade %>% filter(str_detect(Name, "qval")) %>%
  dplyr::rename(qval = val) %>%
  mutate(Name = gsub("_qval", "", .$Name)) %>%
  mutate(NumGenes = str_extract(.$Name, "top[0-9]+")) %>%
  mutate(NumGenes = as.numeric(gsub("top", "", .$NumGenes))) %>%
  mutate(Condition = "CTL") %>%
  mutate(Grade = "low")

CTL_sliced_percentOverlap_lowgrade_BH <- CTL_sliced_percentOverlap_lowgrade %>% filter(str_detect(Name, "FDR")) %>%
  dplyr::rename(FDR = val) %>%
  mutate(Name = gsub("_FDR", "", .$Name)) %>%
  mutate(NumGenes = str_extract(.$Name, "top[0-9]+")) %>%
  mutate(NumGenes = as.numeric(gsub("top", "", .$NumGenes))) %>%
  mutate(Condition = "CTL") %>%
  mutate(Grade = "low")

CTL_sliced_percentOverlap_highgrade <- data.frame("Name" = names(CTL_sliced_percentOverlap_highgrade), 
                                                  "val" = unlist(CTL_sliced_percentOverlap_highgrade))

CTL_sliced_percentOverlap_highgrade_qval <- CTL_sliced_percentOverlap_highgrade %>% filter(str_detect(Name, "qval")) %>%
  dplyr::rename(qval = val) %>%
  mutate(NumGenes = str_extract(.$Name, "top[0-9]+")) %>%
  mutate(NumGenes = as.numeric(gsub("top", "", .$NumGenes))) %>%
  mutate(Name = gsub("_top[0-9]+_qval", "", .$Name)) %>%
  mutate(Condition = "CTL") %>%
  mutate(Grade = "high")

CTL_sliced_percentOverlap_highgrade_BH <- CTL_sliced_percentOverlap_highgrade %>% filter(str_detect(Name, "FDR")) %>%
  dplyr::rename(FDR = val) %>%
  mutate(NumGenes = str_extract(.$Name, "top[0-9]+")) %>%
  mutate(NumGenes = as.numeric(gsub("top", "", .$NumGenes))) %>%
  mutate(Name = gsub("_top[0-9]+_FDR", "", .$Name)) %>%
  mutate(Condition = "CTL") %>%
  mutate(Grade = "high")


CTL_sliced_percentOverlap_all <- bind_rows(left_join(CTL_sliced_percentOverlap_lowgrade_qval,
                                                     CTL_sliced_percentOverlap_lowgrade_BH),
                                           left_join(CTL_sliced_percentOverlap_highgrade_qval,
                                                     CTL_sliced_percentOverlap_highgrade_BH)) %>%
  mutate(Name = gsub("_top[0-9]+", "", .$Name))


FNF_sliced_percentOverlap_lowgrade <- list()
FNF_sliced_percentOverlap_highgrade <- list()
for (option in FNF_options$Name){
  # Read in eGenes for set of covariate options and FDR
  eGenes_FDR <- read_csv(paste0("output/qtl/FNF_", option, "_perm1Mb_FDR.txt")) %>%
    filter(FDR <= 0.05)
  n_slices <- ceiling(nrow(eGenes_FDR)/50) 
  
  for (n in seq(1, n_slices, by = 1)){
    n_genes <- n*50
    slice <- eGenes_FDR %>% slice_min(order_by = FDR, n = n_genes) %>%
      pull(gene_name)
    
    # Calculate overlap with low grade
    percOverlap_lowgrade <- percentOverlap(slice, lowgrade_eGenes$gene_name)
    # Calculate overlap with high grade
    percOverlap_highgrade <- percentOverlap(slice, highgrade_eGenes$gene_name)
    FNF_sliced_percentOverlap_lowgrade[[paste0(option, "_top", n_genes, "_FDR")]] <- percOverlap_lowgrade
    FNF_sliced_percentOverlap_highgrade[[paste0(option, "_top", n_genes, "_FDR")]] <- percOverlap_highgrade
  }
  
  eGenes_qval <- read_csv(paste0("output/qtl/FNF_", option, "_perm1Mb_FDR.txt")) %>%
    filter(qval <= 0.05)
  n_slices <- ceiling(nrow(eGenes_qval)/50)
  
  for (n in seq(1, n_slices, by = 1)){
    n_genes <- n*50
    slice <- eGenes_qval %>% slice_min(order_by = qval, n = n_genes) %>%
      pull(gene_name)
    
    # Calculate overlap with low grade
    percOverlap_lowgrade <- percentOverlap(slice, lowgrade_eGenes$gene_name)
    # Calculate overlap with high grade
    percOverlap_highgrade <- percentOverlap(slice, highgrade_eGenes$gene_name)
    FNF_sliced_percentOverlap_lowgrade[[paste0(option, "_top", n_genes, "_qval")]] <- percOverlap_lowgrade
    FNF_sliced_percentOverlap_highgrade[[paste0(option, "_top", n_genes, "_qval")]] <- percOverlap_highgrade
  }
  
}

FNF_sliced_percentOverlap_lowgrade <- data.frame("Name" = names(FNF_sliced_percentOverlap_lowgrade), 
                                                 "val" = unlist(FNF_sliced_percentOverlap_lowgrade))

FNF_sliced_percentOverlap_lowgrade_qval <- FNF_sliced_percentOverlap_lowgrade %>% filter(str_detect(Name, "qval")) %>%
  dplyr::rename(qval = val) %>%
  mutate(Name = gsub("_qval", "", .$Name)) %>%
  mutate(NumGenes = str_extract(.$Name, "top[0-9]+")) %>%
  mutate(NumGenes = as.numeric(gsub("top", "", .$NumGenes))) %>%
  mutate(Condition = "FNF") %>%
  mutate(Grade = "low")

FNF_sliced_percentOverlap_lowgrade_BH <- FNF_sliced_percentOverlap_lowgrade %>% filter(str_detect(Name, "FDR")) %>%
  dplyr::rename(FDR = val) %>%
  mutate(Name = gsub("_FDR", "", .$Name)) %>%
  mutate(NumGenes = str_extract(.$Name, "top[0-9]+")) %>%
  mutate(NumGenes = as.numeric(gsub("top", "", .$NumGenes))) %>%
  mutate(Condition = "FNF") %>%
  mutate(Grade = "low")

FNF_sliced_percentOverlap_highgrade <- data.frame("Name" = names(FNF_sliced_percentOverlap_highgrade), 
                                                  "val" = unlist(FNF_sliced_percentOverlap_highgrade))

FNF_sliced_percentOverlap_highgrade_qval <- FNF_sliced_percentOverlap_highgrade %>% filter(str_detect(Name, "qval")) %>%
  dplyr::rename(qval = val) %>%
  mutate(NumGenes = str_extract(.$Name, "top[0-9]+")) %>%
  mutate(NumGenes = as.numeric(gsub("top", "", .$NumGenes))) %>%
  mutate(Name = gsub("_top[0-9]+_qval", "", .$Name)) %>%
  mutate(Condition = "FNF") %>%
  mutate(Grade = "high")

FNF_sliced_percentOverlap_highgrade_BH <- FNF_sliced_percentOverlap_highgrade %>% filter(str_detect(Name, "FDR")) %>%
  dplyr::rename(FDR = val) %>%
  mutate(NumGenes = str_extract(.$Name, "top[0-9]+")) %>%
  mutate(NumGenes = as.numeric(gsub("top", "", .$NumGenes))) %>%
  mutate(Name = gsub("_top[0-9]+_FDR", "", .$Name)) %>%
  mutate(Condition = "FNF") %>%
  mutate(Grade = "high")


FNF_sliced_percentOverlap_all <- bind_rows(left_join(FNF_sliced_percentOverlap_lowgrade_qval,
                                                     FNF_sliced_percentOverlap_lowgrade_BH),
                                           left_join(FNF_sliced_percentOverlap_highgrade_qval,
                                                     FNF_sliced_percentOverlap_highgrade_BH)) %>%
  mutate(Name = gsub("_top[0-9]+", "", .$Name))

# Joining and plotting p-value ranked percent overlaps -----------------------------------
sliced_percentOverlap_all <- bind_rows(CTL_sliced_percentOverlap_all, FNF_sliced_percentOverlap_all) %>%
  mutate(Grade = factor(Grade, levels = c("low", "high"))) %>%
  mutate(Name = factor(Name))

ggplot(data = sliced_percentOverlap_all, aes(x = NumGenes, y = FDR, group = Name, color = Name)) +
  geom_line() +
  facet_wrap(~Condition + Grade,
             scales = "free_y")+
  xlab("Top Genes") +
  ylab("Percent Overlap") +
  labs(color = "Legend") +
  ylim(0, 0.8) +
  theme_minimal() + 
  theme(panel.grid.major.y = element_blank(),
        legend.title = element_blank())

ggsave(filename = "output/plots/Steinberg_topGene_percOverlap_qvalue.pdf")

ggsave(filename = "output/plots/Steinberg_topGene_percOverlap_FDR.pdf") 
