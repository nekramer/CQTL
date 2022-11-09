library(tidyverse)
library(RColorBrewer)

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


CTL <- allOptions[[1]] %>% arrange(FDR)
CTL$Name <- factor(CTL$Name, levels = CTL$Name)
FNF <- allOptions[[2]]
FNF$Name <- factor(FNF$Name, levels = CTL$Name)

allOptions <- bind_rows(list(CTL, FNF)) %>% filter(Nk %in% common_Nk)





ggplot(data = allOptions) +
  geom_segment(aes(x = FDR, xend = qval, y = Name, yend = Name), lwd = 0.5) +
  geom_point(aes(x = FDR, y = Name,color = brewer.pal(n = 3, "Paired")[1]), size = 2) +
  geom_point(aes(x = qval, y = Name,color = brewer.pal(n = 3, "Paired")[2]), size = 2) +
  facet_wrap(~Condition,
             scales = "free_y")+
  scale_color_identity() +
  xlab("# eGenes") +
  labs(color = "Legend") +
  theme_minimal() + 
  theme(axis.title.y=element_blank(),
        panel.grid.major.y = element_blank(),
        legend.title = element_blank()) +
  scale_x_continuous(limits = c(1250, 2250), 
                     breaks = c(1250, 1500, 1750, 2000, 2250),
                     labels = c("1250", "1500", "1750", "2000", "2250")) +
  scale_color_identity(guide = "legend",
                       labels = c("q-value", "BH"))
ggsave(filename = "output/plots/covariate_eGenes.pdf")  

