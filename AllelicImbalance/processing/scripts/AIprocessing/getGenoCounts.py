import allel
import sys
import numpy.ma as ma
import numpy as np
import pandas as pd

# sys.argv[1] = path to vcf file (make sure it has its corresponding index file)
# sys.argv[2] = donor name (from RNA)
# sys.argv[3] = path to text file of matching donor names in vcf file
#sys.argv[4] = path to common variants list to subset for

## Grab the donor name in the vcf file
donorNames = pd.read_csv(sys.argv[3], sep = " ")
VCFdonor = donorNames.loc[donorNames.rna == sys.argv[2]].vcf.values[0]

## Read in variant ID and genotype fields for specified donor
callset = allel.read_vcf(sys.argv[1], samples = [VCFdonor], fields = ['ID','GT'])
# Grab genotypes and put into GenotypeArray
gt = allel.GenotypeArray(callset['calldata/GT'])
# Grab variant IDs
variants = callset['variants/ID']

## Grab booleans for which genos are hom and which are het
homs = gt.is_hom()
hets = gt.is_het()

## Used masked arrays to get variant IDs from hom and het geno boolean arrays
# Mask with hets, get homs
homVariants = ma.masked_array(variants, mask = hets)
homVariants = homVariants[~homVariants.mask].data
# Mask with homs, get hets
hetVariants = ma.masked_array(variants, mask = homs)
hetVariants = hetVariants[~hetVariants.mask].data

## Initialize same sized numpy arrays of weights for homVariants (1e-6) and hetVariants (1)
homWeights = np.full_like(homVariants, 1e-6)
hetWeights = np.full_like(hetVariants, 1)

## Combine homVariants/homWeights and hetVariants/hetWeights into dataframes
homDF = pd.DataFrame({'variantID': homVariants, 'genoweight': homWeights})
hetDF = pd.DataFrame({'variantID': hetVariants, 'genoweight': hetWeights})

## Concatenate homDF and hetDF and sort by variantID
allDF = pd.concat([homDF, hetDF], axis = 0)
allDF = allDF.sort_values(by = 'variantID')

## Write to file
# File name
outputName = 'output/AI/' + str(sys.argv[2]) + '/' + str(sys.argv[2]) + '_genoCounts.csv'
allDF.to_csv(outputName, index = False)
