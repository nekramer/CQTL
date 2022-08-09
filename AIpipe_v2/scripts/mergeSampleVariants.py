import pandas as pd
import numpy as np
import sys
import os
import re

# sys.argv[1] = alleleCount file
# sys.argv[2] = list of all variants
# sys.argv[3] = min total number of allele reads both alleles must have to be considered heterozygous
# sys.argv[4] = min number of allele reads each allele must have to be considered heterozygous

# Read in sample's alleleCounts data (this is tab-delimited) 
# Select only "variantID", "refCount", and "altCount" columns
sampleAlleleCounts = pd.read_csv(sys.argv[1], sep = '\t', usecols = ['variantID', 'refCount', 'altCount'])

# Read in unique set of all variants
variants = pd.read_csv(sys.argv[2])

minTotal = int(sys.argv[3])
minAllele = int(sys.argv[4])

# Perform a join with variants to include variants that sample doesn't have
# Joined on variants to preserve variantID order of variants
sampleAlleleCounts_joined = variants.merge(sampleAlleleCounts, on = 'variantID', how = 'outer')
sampleAlleleCounts_joined.set_index('variantID')

# Go through each variant and record a weight value
for variant in sampleAlleleCounts_joined.itertuples():
    if np.isnan(variant.refCount):
        #variant.refCount = 0
        #variant.altCount = 0
        print(sampleAlleleCounts_joined.at[variant.Index, 'refCount'])
        sampleAlleleCounts_joined.at[variant.Index, 'refCount'] = 0
        sampleAlleleCounts_joined.at[variant.Index, 'altCount'] = 0
        sampleAlleleCounts_joined.at[variant.Index, 'weight'] = 0
    elif (variant.refCount + variant.altCount >= minTotal) and (min(variant.refCount, variant.altCount) >= minAllele):
        # Considered a het, so give weight of 1
        sampleAlleleCounts_joined.at[variant.Index, 'weight'] = 1
    else:
        # Considered a hom, so give weight of 1e-6
         sampleAlleleCounts_joined.at[variant.Index, 'weight'] = 1e-6

# Parse sample name from file
sampleName = re.sub("_alleleCounts.csv", "", os.path.basename(sys.argv[1]))
sampleName_split = sampleName.split("_")
donor = sampleName_split[1]
condition = sampleName_split[3]

# Add donor and condition column
sampleAlleleCounts_joined['donor'] = donor
sampleAlleleCounts_joined['condition'] = condition

# Write to file using sample name
outputName = 'output/' + sampleName + '/alleleCounts/' + sampleName + '_alleleCounts_joined.csv'
sampleAlleleCounts_joined.to_csv(outputName, index = False)