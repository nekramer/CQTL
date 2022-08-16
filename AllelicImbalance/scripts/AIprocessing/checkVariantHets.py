import pandas as pd
import sys
import numpy as np
from dask import dataframe as dd

# sys.argv[1]: csv file of concatenated and checked allele counts for all variants

# Function to count how many donor hets a variant has
def variantHets(variantDataGroup):
    # Get one representative row per donor
    variantDataDonors = variantDataGroup.drop_duplicates(['donor'])
    
    # Count number of hets
    numHets = np.sum(variantDataDonors['weight'] == 1)

    return numHets

alleleCounts = dd.read_csv(sys.argv[1])

# Group alleleCounts by variant and count number of het donors per variant
numVariantHets = alleleCounts.groupby('variantID').apply(variantHets, meta = ('x', 'f8'))

# Write Series to csv
numVariantHets.compute().to_csv("output/AI/numVariantHets.csv")