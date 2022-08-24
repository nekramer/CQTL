import pandas as pd
from utils import getRsids
import sys

# sys.argv[1]: LD buddies text file
# sys.argv[2]: OA subtype


# Read in LD buddies dataset
ld_data = pd.read_csv(sys.argv[1], sep = "\t")
# Just grab columsn with seqnames/start
ld_data_positions = ld_data.iloc[:, 1:3]
# Add "chr" prefix to seqnames column to match dbSNP
ld_data_positions['seqnames']= ld_data_positions['seqnames'].apply(lambda x: f"chr{x}")

# Convert to RSIDS using GRCh37 dbSNP vcf 
rsids = getRsids(ld_data_positions, dbSNP = '/proj/phanstiel_lab/References/genomes/GENCODE.GRCh37.p13/dbSNP/dbSNP155.GRCh37.p13.vcf.gz')
ld_data_all = pd.concat([ld_data_positions, rsids], axis = 1)
ld_data_all.to_csv('data/Boer_' + sys.argv[2] + '_LD_rsids.csv', index = False)


