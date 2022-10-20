import pandas as pd
import sys
import numpy as np
from kneed import KneeLocator

# sys.argv[1] = PC file with sdev in last row

def main():
    # Read in 
    PC = pd.read_csv(sys.argv[1])

    # Grab sdev
    sdev = PC.loc[PC["Donor"] == 'sdev'].to_numpy()

    # Get rid of 'sdev' in first index
    sdev = np.delete(sdev, 0)

    # Calculate variance explained by each PC
    varExplained = np.square(sdev)/sum(np.square(sdev))

    # Generate vector for number of PCs
    x = range(1, len(PC.columns))

    # Calculate knee
    kneedle = KneeLocator(x, varExplained, curve = "convex", direction = "decreasing")

    numPCs = kneedle.knee
    return numPCs


