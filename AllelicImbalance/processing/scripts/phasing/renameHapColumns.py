import pandas as pd
import numpy as np
import sys

# sys.argv[1]: path to phased haplotype file
# sys.argv[2]: path to sample file
# sys.argv[3]: chromosome
# sys.argv[4]: vcf prefix

haps = pd.read_csv(sys.argv[1], header = None)

# Read in sample file, skipping first row
samples = pd.read_csv(sys.argv[2], sep = " ", skiprows = [1], usecols = ['ID_1', 'ID_2'])

# Split the strings in each column and just keep the 3rd column (i.e. donor)
col1 = samples.ID_1.str.split("_", expand = True)[3]
col2 = samples.ID_2.str.split("_", expand = True)[3]

# Add '_hap1' and '_hap2' after donor names
col1 = np.core.defchararray.add(col1.to_numpy().astype('U'), '_hap1')
col2 = np.core.defchararray.add(col2.to_numpy().astype('U'), '_hap2')

donors = pd.DataFrame({1: col1, 2: col2})

# Flatten row-wise to get column order
sampleColumns = donors.to_numpy().flatten()

# Initialize first few column names 
firstColumns = np.array(['chrom', 'variantID', 'pos', 'hap1', 'hap2'])

# Concatenate all columns into one array
allColumns = np.concatenate((firstColumns, sampleColumns))

# Set column names
haps.columns = allColumns

haps.to_csv("output/vcf/" + sys.argv[4] + "_nodups_biallelic_chr" + sys.argv[3] + ".phased.allele.named.haps", index = False)