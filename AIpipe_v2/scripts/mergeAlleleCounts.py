import chunk
import pandas as pd
import sys
import os
import re

# sys.argv[1] = list of all variants
# sys.argv[2]...sys.argv[n] = alleleCount file from each sample

# Read in list of variants
variants = pd.read_csv(sys.argv[1])

# Initialize dictionary where variants are the keys
alleleCountsDictionary = dict.fromkeys(variants.variantID, {})

for variant in variants.variantID:
    for file in range(2, len(sys.argv)):
        sample = re.sub("_alleleCounts.csv", "", os.path.basename(sys.argv[file]))
        for line in pd.read_csv(sys.argv[file], sep = "\t", chunksize = 1):
            globals()[sample] = {'ref': None, 'alt': None, 'weight': None}
            if line.variantID.values == variant:
                globals()[sample]['ref'] = line.refCount.values
                globals()[sample]['alt'] = line.altCount.values
                globals()[sample]['weight'] = 1
            else:
                globals()[sample]['ref'] = 0
                globals()[sample]['alt'] = 0
                globals()[sample]['weight'] = 0
            # Add sample dictionary to variant
            alleleCountsDictionary[variant][sample] = globals()[sample]


print(alleleCountsDictionary[variant])
# # Go through all sample files
# for file in range(2, len(sys.argv)):
#     # Get sample name of file
#     sample = re.sub("_alleleCounts.csv", "", os.path.basename(sys.argv[file]))

#     # Go through file line by line
#     for line in pd.read_csv(sys.argv[file], sep = "\t", chunksize = 1):
#         if line.variantID.values in variants.variantID.values:
#             ref = line.refCount.values
#             alt = line.altCount.values
#             weight = 1
#         else:
#             ref = 0
#             alt = 0
#             weight = 0
        
#         alleleCountsDictionary[]
        # If variant in that line is in variants, create sub dict with ref/alt allele counts
        # else, create sub dict with None vlaues 




    

