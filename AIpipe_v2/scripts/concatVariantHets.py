import pandas as pd
import sys

# sys.argv[1]: minHets (needed for file name)
# sys.argv[2]...sys.argv[n]: filtered alleleCounts split files

# Add names of files to list
file_list = []

for i in range(2, len(sys.argv)):
    file_list.append(sys.argv[i])

outputName = 'output/AI/alleleCounts_' + str(sys.argv[1]) + 'hets.csv'
with open(outputName, "a+") as output:
    for file in file_list:
        with open(file, "r") as f:
            for line in f:
                output.write(line)