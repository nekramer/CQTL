import sys
import pandas as pd
import difflib

# Read in samples to remove
removeSamples = pd.read_csv(sys.argv[1], sep = " ")

removePLINK = []
for row in removeSamples.itertuples():
    # Grab fam file for batch
    famFile = pd.read_csv('output/binary/' + row.Batch + '.fam', sep = " ", header = None)

    # Match donor to 2nd column
    donorID = difflib.get_close_matches(row.Donor, famFile.iloc[:,1], 1, cutoff = 0.4)[0]

    # Select FID and within-family id columns
    colIDs = famFile.loc[(famFile[1] == donorID)][[0, 1]]

    removePLINK.append(colIDs)

removePLINK_df = pd.concat(removePLINK)

removePLINK_df.to_csv("remove.plink", sep = " ", header = False, index = False)