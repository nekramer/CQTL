import pandas as pd
import sys
import difflib

# sys.argv[1] = donors.txt file with VCF donors
# sys.argv[2] = string of donors from RNA samplesheet

# Read in file containing sample names from vcf file and use as array
donorFile = pd.read_csv(sys.argv[1], header = None)
vcfDonors = donorFile.iloc[:,0].values

# Convert RNA donor string to array
rnaDonors = sys.argv[2].split(",")

donorOrder = []
# Iterate through rnaDonors and find closest match in vcfDonors
for rnaDonor in rnaDonors:
    vcfDonor = difflib.get_close_matches(rnaDonor, vcfDonors, 1, cutoff = 0.4)
    
    # Append to donorOrder
    donorOrder.append(vcfDonor[0])

# Initialize new donorFile dataframe with vcf and rna
newDonorFile = pd.DataFrame()
newDonorFile['vcf'] = donorOrder
newDonorFile['rna'] = rnaDonors

# Write to file
newDonorFile.to_csv("donors.txt", sep = " ", index = False)