library(googlesheets4)
library(dplyr)

GOOGLE_SHEET <- Sys.getenv("GOOGLE_SHEET")
GSHEET_CLIENT_EMAIL <- Sys.getenv("GSHEET_CLIENT_EMAIL")
GSHEET_PRIVATE_KEY <- Sys.getenv("GSHEET_PRIVATE_KEY")

gs4_auth(email = GSHEET_CLIENT_EMAIL,
        token = GSHEET_PRIVATE_KEY)


donors <- read_sheet(ss = GOOGLE_SHEET,
                     sheet = "Donors")
dna <- read_sheet(ss = GOOGLE_SHEET,
                     sheet = "DNAExtractionsLibraries")
rna <- read_sheet(ss = GOOGLE_SHEET,
                     sheet = "RNA ExtractionsLibraries")

# Donors we have all the data back for
rnaDonor <- rna %>% filter(!is.na(Sequencing_Directory)) %>% pull(Donor) %>% unique()
donors %>% filter(Donor %in% rnaDonor) %>% group_by(Sex) %>% count()

