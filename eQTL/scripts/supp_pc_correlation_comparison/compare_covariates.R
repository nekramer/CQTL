#!/usr/bin/R
suppressPackageStartupMessages(library("argparse"))
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("RColorBrewer"))



# Argument parsing --------------------------------------------------------
parser <- ArgumentParser()
parser$add_argument('--covariateOptions', 
                    help = 'List of covariate options to include in comparison.')
parser$add_argument('--output', default = "covariate_eGenes.pdf",
                    help = 'Output file path and name.')

args <- parser$parse_args()

covariateOptions <- unlist(str_split(args$covariateOptions, " "))

# Get number of significant eGenes ----------------------------------------

sig_eGenes <- list()

for (option in covariateOptions){
  for (condition in c("CTL", "FNF")){
    
    # Get files for all PEER factors
    resultFiles <- list.files(paste0("/pine/scr/n/e/nekramer/iteratePEER_", option, "/output/qtl/"),
               pattern = paste0("^", condition, "_PEER_k.*_genoPC_",
                                option, "_perm1Mb_FDR.txt"))
    
    for (file in resultFiles){
      covar_eGenes <- read_csv(paste0("/pine/scr/n/e/nekramer/iteratePEER_",
                                      option, "/output/qtl/", file))
      
      qval_sig <- covar_eGenes %>% dplyr::filter(qval <= 0.05) %>% nrow()
      FDR_sig <- covar_eGenes %>% dplyr::filter(FDR <= 0.05) %>% nrow()
      
      Nk <- as.numeric(gsub("k", "", unlist(str_split(file, "_"))[3]))
      
      
      if (condition == "CTL"){
        cond = "Control"
      } else if (condition == "FNF"){
        cond = "FN-f"
      }
      
      sig_eGenes[[paste0(condition, "_PEER_Nk", Nk, "_", option)]] <- 
        data.frame("qval" = qval_sig, "FDR" = FDR_sig, "Condition" = cond,
                   "Name" = paste0("PEER_Nk", Nk, "_", option))
    }
    
  }
}


# Order based on FDR and relevel ------------------------------------------

covariate_data <- bind_rows(sig_eGenes) %>% 
  group_by(Condition) %>%
  group_split()


CTL <- covariate_data[[1]] %>% arrange(FDR)
CTL$Name <- factor(CTL$Name, levels = CTL$Name)
FNF <- covariate_data[[2]]
FNF$Name <- factor(FNF$Name, levels = CTL$Name)

covariate_data <- bind_rows(list(CTL, FNF))


# Plot --------------------------------------------------------------------
ggplot(data = covariate_data) +
  geom_segment(aes(x = FDR, xend = qval, y = Name, yend = Name), lwd = 0.5) +
  geom_point(aes(x = FDR, y = Name,color = brewer.pal(n = 3, "Paired")[1]), size = 2) +
  geom_point(aes(x = qval, y = Name,color = brewer.pal(n = 3, "Paired")[2]), size = 2) +
  facet_wrap(~Condition,
             scales = "free_x")+
  scale_color_identity() +
  xlab("# eGenes") +
  labs(color = "Legend") +
  theme_minimal() + 
  theme(axis.title.y=element_blank(),
        panel.grid.major.y = element_blank(),
        legend.title = element_blank()) +
  scale_color_identity(guide = "legend",
                       labels = c("q-value", "BH"))


ggsave(filename = args$output, width = 11, height = 7, units = "in")  