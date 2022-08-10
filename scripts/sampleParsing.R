library(googlesheets4)
library(dplyr)

donors <- read_sheet(ss = "https://docs.google.com/spreadsheets/d/1JwLw9D6rMqhHC9BPrZebAN40Wojo-CqbMdiIPXAzkLo/edit#gid=1699779981",
                     sheet = "Donors")
dna <- read_sheet(ss = "https://docs.google.com/spreadsheets/d/1JwLw9D6rMqhHC9BPrZebAN40Wojo-CqbMdiIPXAzkLo/edit#gid=1699779981",
                     sheet = "DNAExtractionsLibraries")
rna <- read_sheet(ss = "https://docs.google.com/spreadsheets/d/1JwLw9D6rMqhHC9BPrZebAN40Wojo-CqbMdiIPXAzkLo/edit#gid=1699779981",
                     sheet = "RNA ExtractionsLibraries")

# Donors we have all the data back for
rnaDonor <- rna %>% filter(!is.na(Sequencing_Directory)) %>% pull(Donor) %>% unique()
donors %>% filter(Donor %in% rnaDonor) %>% group_by(Sex) %>% count()

