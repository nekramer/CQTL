#!/bin/bash

# This script requires an installation of HOMER with added path to a bash profile for loading.

geneFile=$1
bgFile=$2
outDir=$3

mkdir -p ${outDir}

sbatch -t 72:00:00 --mem=4G --wrap="findMotifs.pl ${geneFile}  human ${outDir} -bg ${bgFile}"