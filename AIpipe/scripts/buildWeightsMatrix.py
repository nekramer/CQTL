import pandas as pd
import sys
from joblib import Parallel, delayed
import multiprocessing

# sys.argv[1] = path to filtered alleleCounts matrix
# sys.argv[2] = path to filtered RNAhets
# sys.argv[3] = weight for homozygotes

def assignColumn(columnName, variants, donordfs, homWeight):

    def assignVariantWeight(variant, rnaHets, homWeight):

        # grab variant het value from rnaHets for donor
        hetValue = rnaHets.loc[rnaHets['variant'] == variant]['rnahet'].item()

        if hetValue == True:
            weight = 1
        else:
            weight = homWeight
        return weight

    # Split column name string for donor
    donor = columnName.split("_")[0]

    # Grab donor set of variants from dictionary
    rnaHets = donordfs[donor]

    num_cores = multiprocessing.cpu_count()
    weights = pd.DataFrame(Parallel(n_jobs = num_cores)(delayed(assignVariantWeight)(variant, rnaHets, homWeight) for variant in variants))
    return(weights)

# Read in alleleCounts matrix to get variants/samples
alleleCounts = pd.read_csv(sys.argv[1], sep = " ")
variants = alleleCounts.index.values

# Read in RNAhets to determine which weights to assign (sort index by variant)
rnaHets = pd.read_csv(sys.argv[2], sep = " ")
# Split into dictionary based on donor (each donor has separate df of variants)
donordfs = dict(tuple(rnaHets.groupby('donor')))

homWeight = float(sys.argv[3])

num_cores = multiprocessing.cpu_count()
weights = Parallel(n_jobs = num_cores)(delayed(assignColumn)(columnName, variants, donordfs, homWeight) for columnName in alleleCounts.columns)

weightMatrix = pd.concat(weights, axis = 1)
weightMatrix.columns = alleleCounts.columns
weightMatrix.index = alleleCounts.index

weightMatrix.to_csv("output/AI/weightMatrix.txt")