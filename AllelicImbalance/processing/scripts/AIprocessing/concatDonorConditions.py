import pandas as pd
import sys

# sys.argv[1] : donor name
# sys.argv[2]...sys.argv[n] = alleleCount files for a donor


# Read in each file and add to list
data_list = []

for i in range(2, len(sys.argv)):
    data = pd.read_csv(sys.argv[i])
    data_list.append(data)

# Vertially concatenate all data into one DataFrame
donorData = pd.concat(data_list, axis = 0)

# File name
outputName = 'output/AI/' + str(sys.argv[1]) + '/' + str(sys.argv[1]) + '_alleleCounts_joined.csv'

donorData.to_csv(outputName, index = False)