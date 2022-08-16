import pandas as pd
import sys

# sys.argv[1]...sys.argv[n] = checked alleleCount file from each sample

# Add names of files to list
file_list = []

for i in range(1, len(sys.argv)):
    file_list.append(sys.argv[i])

with open('output/AI/alleleCounts.csv', "a+") as output:
    firstFile = True
    for file in file_list:
        with open(file, "r") as f:
            if not firstFile:
                next(f)
            for line in f:
                output.write(line)
        firstFile = False
