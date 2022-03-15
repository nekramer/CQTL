import pandas as pd
import sys, csv

## Read in data file
data = pd.read_csv(sys.argv[1], sep = "\t", header = None)

setA1_join = []
## Read in setA1 file line by line, and if it matches a position in data, write to dataframe
with open(sys.argv[2]) as setA1:
    for line in setA1:
        pos = line.split()[0]
        chrom = int(pos.split(":")[0])
        chromData = data.loc[data[0] == chrom]
        bp = int(pos.split(":")[1])
        bpData = chromData.loc[chromData[3] == bp]
        if not bpData.empty:
            setA1_join.append([bpData[1].to_string(index = False).strip(), line.split()[1]])

setA1_join = pd.DataFrame(setA1_join)
## Drop duplicate variant ids
setA1_join = setA1_join.drop_duplicates(0)

setA1_join.to_csv(sys.argv[3], header = False, index = False, sep = " ")
