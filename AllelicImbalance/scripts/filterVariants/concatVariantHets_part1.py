import pandas as pd
import sys

# sys.argv[1]: minHets (needed for file name)
# sys.argv[2]: subset group name (needed for file name)
# sys.argv[3]...sys.argv[n]: filtered alleleCounts split files

# Add names of files to list
file_list = []

for i in range(3, len(sys.argv)):
    file_list.append(sys.argv[i])

outputName = 'output/AI/alleleCounts_' + str(sys.argv[1]) + 'hets' + str(sys.argv[2]) + '.csv'
with open(outputName, "a+") as output:
    # Write header
    output.write('variantID,refCount,altCount,donor,condition,weight\n')

    for file in file_list:
        with open(file, "r") as f:
            for line in f:
                output.write(line)