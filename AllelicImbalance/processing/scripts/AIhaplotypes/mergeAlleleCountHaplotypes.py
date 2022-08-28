import pandas as pd
import numpy as np
import sys

# sys.argv[1]: path to subsetted allele counts data
# sys.argv[2]: path to phased haplotype data
# sys.argv[3]: chromosome


def getDonorInfo(donorGroup):
    
    # Removing donor AM7352 since we don't have phasing info on this one
    if donorGroup.name == 'AM7352':
        return(pd.DataFrame())

    else: 
        # Just grab haplotype columns for that donor
        subset = donorGroup[['chrom', 'pos', 'variantID', 'refAllele', 'altAllele', 'refCount', 'altCount', 'donor', 'condition', str(donorGroup.name) + '_hap1', str(donorGroup.name) + '_hap2']]

        # Filter snps that are homozygous, first based on haplotypes, then based on counts
        hetsOnly = subset[subset[str(donorGroup.name) + '_hap1'] != subset[str(donorGroup.name) + '_hap2']]
        hetsOnly = hetsOnly[(hetsOnly['refCount'] >= 2) & (hetsOnly['altCount'] >= 2) & (hetsOnly['refCount'] + hetsOnly['altCount'] >= 10)]

        # Make hap1 as refCount and hap2 as altCount
    
            # If refAllele does not equal hap1, we'll get a 1; if refAllele is the same as hap1, we'll get a 0
        countIndices = (hetsOnly['refAllele'] != hetsOnly[str(donorGroup.name) + '_hap1']).astype(int).to_numpy()
            # Use these ints to grab either refCount or altCount (not equal = 1, grab altCount for hap1;
            # equal = 0, grab refCount for hap1)
        countIndices = np.reshape(countIndices, (-1, 1)) # Using 2D numpy arrays to do this
        hap1Counts = np.take_along_axis(hetsOnly[['refCount', 'altCount']].to_numpy(), countIndices, axis = 1)
    
        # Get total counts for ASEP
        totalCounts = hetsOnly['refCount'] + hetsOnly['altCount']

        # Assembly new DataFrame with required info for ASEP
        ASEP_df = pd.DataFrame({'id':  hetsOnly['donor'].values, 'group': hetsOnly['condition'].values, 'snp': hetsOnly['variantID'].values, 'ref': hap1Counts.flatten(), 'total': totalCounts})

        return(ASEP_df)


alleleCounts = pd.read_csv(sys.argv[1])
haps = pd.read_csv(sys.argv[2])

# Group by donor and assign variant ref counts based on hap1 across donors
ASEP_df_allDonors = alleleCounts.merge(haps, on = ['chrom', 'pos']).groupby('donor').apply(getDonorInfo)

ASEP_df_allDonors.to_csv("output/AI/chr" + sys.argv[3] + "_ASEP.csv", index = False)
