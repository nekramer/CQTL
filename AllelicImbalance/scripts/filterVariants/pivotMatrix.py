import pandas as pd
import sys

# sys.argv[1]: subsetted allele counts data (based on a number of hets)

# Read in allele counts
alleleCounts = pd.read_csv(sys.argv[1])
alleleCounts = alleleCounts.rename(columns = {'refCount': 'ref', 'altCount': 'alt'})

## Counts matrix
countMatrix = alleleCounts.pivot(columns = ['donor','condition'], values = ['ref','alt'], index = 'variantID')

## Weights matrix
# Copy weight into refWeight and altWeight for appropriate pivoting
alleleCounts['altWeight'] = alleleCounts.loc[:,'weight']
alleleCounts = alleleCounts.rename(columns = {'weight': 'refWeight'})
# Pivot into weights matrix
weightMatrix = alleleCounts.pivot(columns = ['donor','condition'], values = ['refWeight', 'altWeight'], index = 'variantID')

## Put columns into colData
colData = pd.DataFrame.from_records(data = countMatrix.columns, columns = ['allele', 'donor', 'condition'])

## Write to files
countMatrix.to_csv('output/AI/alleleCountsMatrix.csv')
weightMatrix.to_csv('output/AI/weightsMatrix.csv')
colData.to_csv('output/AI/colData.csv', index = False)

