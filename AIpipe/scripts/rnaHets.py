import pandas as pd
import numpy as np
import sys
import re

# sys.argv[1] = combined alleleCounts matrix
# sys.argv[2] = donors
# sys.argv[3] = min total number of allele reads both alleles must have to be considered heterozygous
# sys.argv[4] = min number of allele reads each allele must have to be considered heterozygous

def getHet(variantRow, donor, totalThreshold, hetThreshold, df):

    # Separate name to find options in condition/bioreps slots
    conds = []
    reps = []
    for name in variantRow.index:
        conds.append(name.split("_")[1])
        reps.append(name.split("_")[2])

    conds = np.unique(conds)
    reps = np.unique(reps)

    genos = []
    alts = []
    zero = False
    for cond in conds:
        for rep in reps:
            name = re.compile(donor + "_" + cond + "_" + rep + ".")
            if len(list(filter(name.match, variantRow.index))) > 0:
                columns = list(filter(name.match, variantRow.index))
                pair = variantRow[columns]
                # Check if any allele counts are 0
                if zero == False:
                    if any(pair == 0):
                        zero = True
                # Get alt allele counts
                alt = re.compile(donor + "_" + cond + "_" + rep + "_alt")
                altCounts = pair[list(filter(alt.match, pair.index))]
                alts.append(altCounts.to_numpy()[0])
                if sum(pair) >= totalThreshold and min(pair) >= hetThreshold:
                    genos.append(True)
                else:
                    genos.append(False)

    if all(genos):
        het = True
    else:
        het = False

    df.append([variantRow.name, donor, het, max(alts), zero])
    return None

# Read in alleleCount data
alleleCounts = pd.read_csv(sys.argv[1], sep = " ")

# Convert donor/condition strings to vectors
donors = sys.argv[2].split(",")

# Get column names
columnNames = alleleCounts.columns

totalThreshold = int(sys.argv[3])
hetThreshold = int(sys.argv[4])

# Initialize a list for compiling variant/donor/het info
hetInfo = []

for donor in donors:
    # Grep for columns for that donor (all conditions/bioreps/alleles)
    name = re.compile(donor + "_" + ".")
    columns = list(filter(name.match, columnNames))

    # Subset alleleCounts for those columns
    donor_alleleCounts = alleleCounts[columns]

    # Apply getHet function to every variant row
    donor_alleleCounts.apply(lambda row: getHet(row, donor, totalThreshold, hetThreshold, hetInfo), axis = 1)
     

# Concat hetInfo into one pandas dataframe
hetInfo = pd.DataFrame(hetInfo, columns = ["variant", "donor", "rnahet", "maxAlt", "anyzero"])

# Write to file
hetInfo.to_csv("output/AI/RNAhets.csv", index = False)

