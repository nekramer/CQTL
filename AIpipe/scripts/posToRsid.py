from pysam import VariantFile
import pandas as pd
import sys

# sys.argv[1] = path to SnpDb vcf file (make sure it has its corresponding index file)
# sys.argv[2] = path to file of variant locations (chr in first column, pos in second column)
# sys.argv[3] = output

# Read in snpDb vcf file
snpInfo = VariantFile(sys.argv[1])

# Read in variant file locations
variants = pd.read_csv(sys.argv[2], sep = " ")

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


pd.DataFrame(rsids).to_csv(sys.argv[3], header = ["rsid"], index = False)