library(googlesheets4)
library(dplyr)

GOOGLE_SHEET <- "*******************"

donors <- read_sheet(ss = GOOGLE_SHEET,
                     sheet = "Donors")
dna <- read_sheet(ss = GOOGLE_SHEET,
                     sheet = "DNAExtractionsLibraries")
rna <- read_sheet(ss = GOOGLE_SHEET,
                     sheet = "RNA ExtractionsLibraries")

# Donors we have all the data back for
rnaDonor <- rna %>% filter(!is.na(Sequencing_Directory)) %>% pull(Donor) %>% unique()
donors %>% filter(Donor %in% rnaDonor) %>% group_by(Sex) %>% count()

# To be extracted
extractedRNA <-  rna %>% filter(is.na(Sequencing_Directory)) %>% pull(Donor) %>% unique()
donors %>% filter(Donor %in% extractedRNA) %>% group_by(Sex) %>% count()