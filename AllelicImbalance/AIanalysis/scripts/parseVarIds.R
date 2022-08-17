source('scripts/utils.R')

# Parsing variantIDs ------------------------------------------------------

## CTL

# Load data
load('data/2022-08-17_AIresCTL.rda')

# Get GRanges of variant positions
AIresCTL_ranges <- varID_to_GRanges(rownames(notNA_resCTL))

write.csv(AIresCTL_ranges, file = 'data/AIresCTL.csv',
          quote = FALSE, row.names = FALSE )


## FNF
load('data/2022-08-17_AIresFNF.rda')

AIresFNF_ranges <- varID_to_GRanges(rownames(notNA_resFNF))

write.csv(AIresFNF_ranges, file = 'data/AIresFNF.csv',
          quote = FALSE, row.names = FALSE)
