import pandas as pd
import sys

# sys.argv[1]: txt file containing names of split concatenated/check allele counts files 
# sys.argv[2]: csv file that quantifies the number of heterozygote donors for each variant
# sys.argv[3]: minimum het threshold for a variant

## Read in number of variant hets
numVariantHets = pd.read_csv(sys.argv[2])

## Filter numVariantHets for the variants that satisfy the minimum het threshold
hetVars = numVariantHets.loc[numVariantHets.x >= int(sys.argv[3])]

outputName = str(sys.argv[1]) + '_' + str(sys.argv[3]) + 'hets'
with open(outputName, "a+") as output:
    with open(sys.argv[1], "r") as f:
        # Check for first split file for header
        if str(sys.argv[1]) == "alleleCounts_splitaa":
            firstLine = True
        else: 
            firstLine = False

        for line in f:
            if firstLine:
                output.write(line)
            else:
                variant = line.split(",")[0]
                if variant in hetVars.variantID.tolist():
                    output.write(line)
            firstLine = False    

