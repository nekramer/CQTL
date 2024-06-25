#!/bin/bash

# This script requires an installation of HOMER with added path to a bash profile for loading.

peakFile=$1
bgFile=$2
outDir=$3

mkdir -p ${outDir}

sbatch -t 72:00:00 --mem=4G --wrap="findMotifsGenome.pl ${peakFile} hg38 ${outDir} -bg ${bgFile}"