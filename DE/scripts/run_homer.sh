#!/bin/bash

# This script requires an installation of HOMER with added path to a bash profile for loading.

geneFile=$1
outDir=$2

mkdir -p ${outDir}

# Grab gene ids (column 6) and write to temp file
awk -F ',' '{print $6}' ${geneFile} > ${outDir}/gene_ids.txt 
# Remove header
sed -i '1d' ${outDir}/gene_ids.txt 


sbatch -t 72:00:00 --mem=4G --wrap="findMotifs.pl ${outDir}/gene_ids.txt  human ${outDir}"