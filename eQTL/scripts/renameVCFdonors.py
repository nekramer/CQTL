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
# Iterate through vcfDonors and find closest match in rnaDonors
for vcfDonor in vcfDonors:
    rnaDonor = difflib.get_close_matches(vcfDonor, rnaDonors, 1, cutoff = 0.4)
    
    # Append to donorOrder
    donorOrder.append(rnaDonor[0])

newHeader = pd.DataFrame(donorOrder)
newHeader.to_csv("samples.txt", sep = " ", index = False, header = False)