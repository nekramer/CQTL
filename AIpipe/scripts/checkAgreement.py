import pandas as pd
import numpy as np
import sys
from joblib import Parallel, delayed
import multiprocessing

# sys.argv[1] = path to genoHets
# sys.argv[2] = path to rnaHets

# Read in geno hets data 
geno = pd.read_csv(sys.argv[1], sep = ",")

# Read in rna hets data
rna = pd.read_csv(sys.argv[2], sep = ",")

# Join geno hets and rna hets data
merge = geno.merge(rna, on = ['variant', 'donor'], how = 'inner')

# Break into a dictionary of dataframes, separated by variant
variantdfs = dict(tuple(merge.groupby('variant')))

def checkMatch(key, variant):
    # Count number of matching

    matching = np.sum(variant.apply(lambda x: x.genohet == x.rnahet, axis = 1))

    # Count number of mismatching
    mismatching = np.sum(variant.apply(lambda x: x.genohet != x.rnahet, axis = 1))

    return matching, mismatching

num_cores = multiprocessing.cpu_count()
stats = Parallel(n_jobs = num_cores)(delayed(checkMatch)(key, variant) for key, variant in variantdfs.items())

matching = sum(i for i, j in stats)
mismatching = sum(j for i, j in stats)
total = matching + mismatching

matching_perc = (matching/total)*100
mismatching_perc = (mismatching/total)*100

print("Matching: " + str(matching) + " (" + str(matching_perc) + "%)")
print("Mismatching: " + str(mismatching) + " (" + str(mismatching_perc) + "%)")