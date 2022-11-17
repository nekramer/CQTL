library(tidyverse)
source("scripts/utils.R")


lowgrade <- read_delim("/proj/phanstiel_lab/External/public/Steinberg_2021/eQTL_LowGrade_FastQTL.txt", 
                              delim = "\t", 
                              col_types = "ccddddddddcdddddddddcc") %>% 
  separate(Gene, into = c("gene_name", "gene_id"), sep = "_")

highgrade <- read_delim("/proj/phanstiel_lab/External/public/Steinberg_2021/eQTL_HighGrade_FastQTL.txt", 
                               delim = "\t", 
                               col_types = "ccddddddddcdddddddddcc") %>%
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
covariateNames <- allOptions[[1]] %>% filter(Nk %in% common_Nk) %>% pull(Name)



eQTLs <- list()
for (condition in conditions){
  for (option in covariateNames){
    
    num_variants <- read_csv(paste0("/pine/scr/n/e/nekramer/CQTL/eQTL/output/qtl/", 
                                condition,"_",
                                option, "_nom1Mb_final.txt")) %>% nrow()
    eQTLs[[paste0(condition, "_", option)]] <- num_variants
    
    
    
    
    
    
    
  }
}

eQTLs_df <- data.frame(num_eQTLs = unlist(eQTLs), Name = names(eQTLs)) %>%
  mutate(Condition = factor(unlist(lapply(str_split(.$Name, "_", n = 2), `[[`, 1)))) %>%
  mutate(Name = str_remove(unlist(lapply(str_split(.$Name, "_", n = 2), `[[`, 2)), "_FDR")) %>%
  group_by(Condition) %>%
  group_split()

eQTLs_CTL <- eQTLs_df[[1]] %>% arrange(num_eQTLs)
eQTLs_CTL$Name <- factor(eQTLs_CTL$Name, levels = eQTLs_CTL$Name)
eQTLs_FNF <- eQTLs_df[[2]]
eQTLs_FNF$Name <- factor(eQTLs_FNF$Name, levels = eQTLs_CTL$Name)

eQTLs_all <- bind_rows(eQTLs_CTL, eQTLs_FNF) 




ggplot(data = eQTLs_all) +
    geom_point(aes(x = num_eQTLs, y = Name), size = 2) +
    facet_wrap(~Condition,
               scales = "free_y")+
    xlab("Number of eQTLs") +
    labs(color = "Legend") +
    theme_minimal() +
    theme(axis.title.y=element_blank(),
          panel.grid.major.y = element_blank(),
          legend.title = element_blank())

ggsave(filename = "output/plots/covariate_eQTLs.pdf")
