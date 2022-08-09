import pandas as pd
import sys
import os
import re

# sys.argv[1] = list of all variants
# sys.argv[2] alleleCount file

# Read in unique set of all variants
variants = pd.read_csv(sys.argv[1])

# Read in sample's alleleCounts data (this is tab-delimited) 
# Select only "variantID", "refCount", and "altCount" columns
sampleAlleleCounts = pd.read_csv(sys.argv[2], sep = '\t', usecols = ['variantID', 'refCount', 'altCount'])

# Perform a join with variants to include variants that sample doesn't have
# Joined on variants to preserve variantID order of variants
sampleAlleleCounts_joined = variants.merge(sampleAlleleCounts, on = 'variantID', how = 'outer')

# Parse sample name from file to write new file
sampleName = re.sub("_alleleCounts.csv", "", os.path.basename(sys.argv[2]))
outputName = 'output/' + sampleName + '/alleleCounts/' + sampleName + '_alleleCounts_joined.csv'

# Write to file
sampleAlleleCounts_joined.to_csv(outputName)