from pysam import VariantFile
import pandas as pd
import sys

# sys.argv[1] = path to vcf file (make sure it has its corresponding index file)
# sys.argv[2] = path to text file of matching donor names in vcf file

# Read in vcf for geno data
geno = VariantFile(sys.argv[1])

# Read in file indicating corresponding vcf/rna donor names
donorNames = pd.read_csv(sys.argv[2], sep = " ")
geno_donors = list(geno.header.samples)

hetInfo = []

# Iterate through all variants
for variant in geno.fetch():
    # Iterate through all donors
    for d in geno_donors:
        # Get corresponding name of donor in rna
        donor = donorNames.loc[donorNames['vcf'] == d]['rna'].values[0]
        # Grab genotype
        genotype = variant.samples[d]['GT']
        # If alleles are the same consider homozygous (weight of 1e-6), if different consider heterozygous (weight of 1)
        if genotype[0] == genotype[1]:
            weight = 1e-6
        else:
            weight = 1
        hetInfo.append([variant.id, donor, weight])

# Concat hetInfo into one pandas dataframe
hetInfo = pd.DataFrame(hetInfo, columns = ["variant", "donor", "weight"])

# Write to file
#hetInfo.to_csv("output/AI/genohets.csv", index = False)
hetInfo.to_pickle('output/AI/genohets.pkl')
    

    
    







