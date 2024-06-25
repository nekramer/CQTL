library(tidyverse)
source("../utils.R")
library(rrvgo)

## Write high confidence response eGenes to file
pbs_highconf <- read_csv("data/reqtl/CTL_sig01_beta_donor_reQTLs_PEER_k20_genoPC.csv") |> 
  dplyr::select(gene_id) |> 
  write_csv(file = "data/reqtl/homer/PBS/PBS_highconf_reGenes.csv", col_names = FALSE)

fnf_highconf <- read_csv("data/reqtl/FNF_sig01_beta_donor_reQTLs_PEER_k22_genoPC.csv") |> 
  dplyr::select(gene_id) |> 
  write_csv(file = "data/reqtl/homer/FNF/FNF_highconf_reGenes.csv", col_names = FALSE)

all_highconf <- bind_rows(pbs_highconf, fnf_highconf) |> 
  distinct() |> 
  write_csv(file = "data/reqtl/homer/ALL/ALL_highconf_reGenes.csv", col_names = FALSE)

# Launch homer
system("scripts/highconf_reGenes_homer.sh data/reqtl/homer/PBS/PBS_highconf_reGenes.csv data/reqtl/homer/PBS")
system("scripts/highconf_reGenes_homer.sh data/reqtl/homer/FNF/FNF_highconf_reGenes.csv data/reqtl/homer/FNF")
system("scripts/highconf_reGenes_homer.sh data/reqtl/homer/ALL/ALL_highconf_reGenes.csv data/reqtl/homer/ALL")

pbs_regene_go <- read_delim("data/reqtl/homer/PBS/biological_process.txt") |> 
  mutate(pval = exp(1)^logP) |> 
  filter(pval < 0.01)

pbs_regene_go <- reduceGO(pbs_regene_go,
                          category = "PBS")

pbs_regene_kegg <- read_delim("data/reqtl/homer/PBS/kegg.txt") |> 
  mutate(pval = exp(1)^logP) |> 
  filter(pval < 0.01) |> 
  distinct(Term, .keep_all = TRUE) |> 
  mutate(`-log10pval` = -log10(pval))



fnf_regene_go <- read_delim("data/reqtl/homer/FNF/biological_process.txt") |> 
  mutate(pval = exp(1)^logP) |> 
  filter(pval < 0.01)

fnf_regene_go <- reduceGO(fnf_regene_go,
                          category = "FN-f")

fnf_regene_kegg <- read_delim("data/reqtl/homer/FNF/kegg.txt") |> 
  mutate(pval = exp(1)^logP) |> 
  filter(pval < 0.01) |> 
  distinct(Term, .keep_all = TRUE) |> 
  mutate(`-log10pval` = -log10(pval))


all_regene_go <- read_delim("data/reqtl/homer/ALL/biological_process.txt") |> 
  mutate(pval = exp(1)^logP) |> 
  filter(pval < 0.01)

all_regene_go <- reduceGO(all_regene_go,
                          category = "ALL")

all_regene_kegg <- read_delim("data/reqtl/homer/ALL/kegg.txt") |> 
  mutate(pval = exp(1)^logP) |> 
  filter(pval < 0.01) |> 
  distinct(Term, .keep_all = TRUE) |> 
  mutate(`-log10pval` = -log10(pval))


## Format and write to table
all_regene_kegg_table <- all_regene_kegg |> 
  dplyr::select(-logP, -pval, -`Entrez Gene IDs`) |> 
  relocate(`-log10pval`, .after = Enrichment) |> 
  arrange(desc(`-log10pval`))

write_csv(all_regene_kegg_table, file = "tables/SupTable10.csv")
ss <- gs4_create(name = "SupTable10")
write_sheet(all_regene_kegg_table,
            ss, sheet = "KEGG pathways - High-confidence PBS-specific and FN-f-response eQTLs")


drive_mv(file = "SupTable10", path = as_dribble("chQTL/chQTL paper/Figures and Tables"))

