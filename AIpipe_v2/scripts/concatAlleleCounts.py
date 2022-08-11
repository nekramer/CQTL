import pandas as pd
import sys

# sys.argv[1]...sys.argv[n] = alleleCount file from each sample

# Read in each file and add to list
data_list = []

for i in range(1, len(sys.argv)):
    data = pd.read_csv(sys.argv[i])
    data_list.append(data)

# Vertially concatenate all data into one DataFrame
allData = pd.concat(data_list, axis = 0)

allData.to_pickle('output/AI/alleleCounts.pkl')