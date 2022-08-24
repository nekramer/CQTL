# utils.py>
import allel
import pandas as pd
# Function to convert variant positions to variant rsids based on a dbSNP reference vcf
def getRsids(variantPositions, dbSNP = '/proj/phanstiel_lab/References/genomes/GENCODE.GRCh38.p13/dbSNP/dbSNP155.GRCh38.p13.vcf.gz'):
    
    def perVariantRSID(variant):
        regionString = str(variant[0]) + ':' + str(variant[1]) + '-' + str(variant[1])
        rsid = allel.read_vcf(dbSNP, region = regionString, fields = ['ID'])
        return(rsid['variants/ID'])
        
    # Make sure input is DataFrame
    variantPositions = pd.DataFrame(variantPositions)

    # Apply perVariantRSID to every row
    variantRSIDs = variantPositions.apply(perVariantRSID, axis = 1)
    return(variantRSIDs)
