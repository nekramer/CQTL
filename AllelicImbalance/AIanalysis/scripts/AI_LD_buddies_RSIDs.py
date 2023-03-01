import pandas as pd
from utils import getRsids
import sys

LD_buddies = sys.argv[1]
snp_name = sys.argv[2]
condition = sys.argv[3]

# Read in output from PLINK and grab SNPB chromosome and position
LD_buddies_positions = pd.read_csv(LD_buddies, delimiter = r"\s+", usecols = ['CHR_B', 'BP_B'])
LD_buddies_positions['CHR_B'] = 'chr' + LD_buddies_positions['CHR_B'].astype(str)


LD_buddies_rsIDs = getRsids(LD_buddies_positions)

LD_buddies_positions_rsIDs = pd.concat([LD_buddies_positions, LD_buddies_rsIDs], axis = 1)
LD_buddies_positions_rsIDs.to_csv('data/' + snp_name + '_' + condition + '_ld_rsids.csv', index = False)
