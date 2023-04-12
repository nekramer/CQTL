from utils import getRsids
import pandas as pd
import sys

# sys.argv[1]: GWAS summary stats
# sys.argv[2]: OA subtype
# sys.argv[3]: chrom


# Read in summary stats
summary_stats = pd.read_csv(sys.argv[1])

# Grab positions
summary_stats_positions = summary_stats.iloc[:, 7:9]


# Convert to RSIDS using GRCh37 dbSNP vcf 
rsids = getRsids(summary_stats_positions, dbSNP = '/proj/phanstiel_lab/References/genomes/GENCODE.GRCh37.p13/dbSNP/dbSNP155.GRCh37.p13.vcf.gz')
summary_stats_all = pd.concat([summary_stats_positions, rsids], axis = 1)
summary_stats_all.to_csv('data/raw/' + sys.argv[2] + '_chr' + str(sys.argv[3]) + 'summary05rsid.csv', index = False)
