import pandas as pd
import sys
import numpy as np
from joblib import Parallel, delayed
import multiprocessing

# sys.argv[1] = path to genoHets
# sys.argv[2] = path to rnaHets
# sys.argv[3] = threshold of number of heterozygotes a variant must have to remain in the analysis

# Read in geno hets data 
geno = pd.read_csv(sys.argv[1], sep = ",")

# Read in rna hets data
rna = pd.read_csv(sys.argv[2], sep = ",")

# Join geno hets and rna hets data
merge = geno.merge(rna, on = ['variant', 'donor'], how = 'inner')

# Break into a dictionary of dataframes, separated by variant
variantdfs = dict(tuple(merge.groupby('variant')))

def checkVariant(key, variant):
    # Check that all genohet rows == rnahet rows
    if all(variant.apply(lambda x: x.genohet == x.rnahet, axis = 1)):

        # Check the number of rnahet == True / anyzero == True donors
        zeroCountHetDonors = np.sum(variant.apply(lambda x : x.rnahet == True and x.anyzero == True, axis = 1))

        # Check the number of rnahet == False / maxAlt >= 10 donors
        altCountHomDonors = np.sum(variant.apply(lambda x : x.rnahet == False and x.maxAlt >= 10, axis = 1))
        if zeroCountHetDonors >= 7 or altCountHomDonors >= 1:
            return key
        else:
            return None
    else:
        return key

num_cores = multiprocessing.cpu_count()
remove_variants = list(filter(None, Parallel(n_jobs = num_cores)(delayed(checkVariant)(key, variant) for key, variant in variantdfs.items())))

# Clear out variant dictionary entries in remove_variants
for var in remove_variants:
    variantdfs.pop(var)

# Now go through to check that each remaining variant has at least 5 heterozygous donors, and if not, add to removal list
def hetThreshold(key, variant, threshold):
    # Count number of het donors
    hets = np.sum(variant.apply(lambda x : x.rnahet == True, axis = 1))
    if hets < threshold:
        return key
    else:
        return None

threshold = int(sys.argv[3])
remove_variants2 = list(filter(None, Parallel(n_jobs = num_cores)(delayed(hetThreshold)(key, variant, threshold) for key, variant in variantdfs.items())))

total_remove = remove_variants + remove_variants2

# Write removal list to file
total_remove = pd.DataFrame(total_remove, columns = ["variant"])
total_remove.to_csv("output/AI/removeVariants.csv", index = False)