import pandas as pd
import sys

# sys.argv[1]: txt file containing names of split concatenated/check allele counts files 
# sys.argv[2]: csv file that quantifies the number of heterozygote donors for each variant
# sys.argv[3]: minimum het threshold for a variant

## Read in split file
if sys.argv[1] == "alleleCounts_splitaa":
    # First file already has header
    splitFile = pd.read_csv(sys.argv[1])
    header = True
else:
    # Other files don't have headers, adding column names
    splitFile = pd.read_csv(sys.argv[1], header = None)
    splitFile.columns = ['variantID', 'refCount', 'altCount', 'donor', 'condition', 'weight']
    header = False


## Read in number of variant hets
numVariantHets = pd.read_csv(sys.argv[2])

## Filter numVariantHets for the variants that satisfy the minimum het threshold
hetVars = numVariantHets.loc[numVariantHets.x >= int(sys.argv[3])]

## Select rows in splitFile whose 'variantID' is in the 'variantID' of hetVars
filteredSplit = splitFile[splitFile['variantID'].isin(hetVars.variantID)]

## Write to file
outputName = str(sys.argv[1]) + '_' + str(sys.argv[3]) + 'hets'
filteredSplit.to_csv(outputName, header = header, index = False)