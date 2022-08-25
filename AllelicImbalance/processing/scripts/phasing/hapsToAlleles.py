import pandas as pd
import sys

# sys.argv[1]: path to phased haplotype file
# sys.argv[2]: chromosome
# sys.argv[3]: vcf prefix

def convertRow(varRow):
    # Replace 0s with first allele
    varRow = varRow.replace(to_replace = 0, value = varRow[3])

    # Replace 1s with second allele
    varRow = varRow.replace(to_replace = 1, value = varRow[4])
    return(varRow)

haps = pd.read_csv(sys.argv[1], sep = " ", header = None)

newHaps = haps.apply(convertRow, axis = 1)

newHaps.to_csv("output/vcf/" + sys.argv[3] + "_nodups_biallelic_chr" + sys.argv[2] + ".phased.allele.haps", index = False, header = None)