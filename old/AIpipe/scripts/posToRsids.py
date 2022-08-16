from pysam import VariantFile
import pandas as pd
import argparse


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('variantLocations', type = argparse.FileType(), 
    help = 'path to space-delimited file of variant locations with 2 columns: "chr" for chromosome name and "pos" for position.')
parser.add_argument('dbSNP',
    help = 'path to dbSNP file to be used for finding variant positions and rsids. The corresponding index must be in the same parent directory.')
parser.add_argument('output', 
    help = 'Desired output file name.')
args = parser.parse_args()

# Read in variant file locations
variants = pd.read_csv(args.variantLocations, sep = " ")

# Read in dbSNP vcf file
snpInfo = VariantFile(args.dbSNP)

rsids = []

for var in variants.itertuples():
    chr = var.chr
    pos = var.pos

    # Use .fetch to grab rsid based on position
    rsid = list(snpInfo.fetch(contig = chr, start = int(pos) - 1, stop = int(pos)))
    if len(rsid) > 0:
        id = rsid[0].info["RS"]
        rsids.append("rs" + str(id))
    else:
        id = chr + ":" + str(pos)
        rsids.append(id)

total = pd.concat([variants, pd.DataFrame(rsids)], axis = 1)
total.to_csv(args.output, header = ["chr", "pos", "rsid"], index = False, sep = " ")
