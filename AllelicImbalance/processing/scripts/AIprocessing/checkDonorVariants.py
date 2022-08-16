import pandas as pd
import sys

# sys.argv[1]: joined alleleCount file from a donor (includes a donor's CTL and FNF)
# sys.argv[2]: quantified genotyping from a donor
# sys.argv[3]: donor

# Function to check whether a donor's CTL, FNF, and geno calls match for a variant
def checkVariant(variantDataGroup):
    # Concatenate weight and genoweight column into one and check how many unique weight values there are
    if pd.concat([variantDataGroup.weight, variantDataGroup.genoweight]).nunique() != 1:
        newweight = 0
    else:
        newweight = variantDataGroup.weight.unique()[0]
    return newweight

# Read in RNA allele counts info and genotyping info from donor
sampleAlleleCounts = pd.read_csv(sys.argv[1])
sampleGeno = pd.read_csv(sys.argv[2])

# Merge data based on variantID
mergedData = sampleAlleleCounts.merge(sampleGeno, how = 'left', on = 'variantID')

# Group merged data by variant and check each variant, returning updated weights
checkedWeights = mergedData.groupby('variantID').apply(checkVariant).rename('newweight')

# Merge new weights with sampleAlleleCounts data
updatedWeightData = sampleAlleleCounts.merge(checkedWeights, right_index = True, left_on = 'variantID')

# Drop old weight column from data and rename newweight column
updatedWeightData.drop(['weight'], axis = 1, inplace = True)
updatedWeightData.rename(columns = {'newweight': 'weight'}, inplace = True)

# Write to file
outputName = 'output/AI/' + str(sys.argv[3]) + '/' + str(sys.argv[3]) + '_alleleCounts_checked.csv'
updatedWeightData.to_csv(outputName, index = False)