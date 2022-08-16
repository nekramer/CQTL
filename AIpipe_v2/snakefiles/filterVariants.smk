#!/usr/bin/env python3

import pandas as pd
import os, re


import numpy as np
import math

## Load config file
configfile: "config/config_filterVariants.yaml"

## Get split file names in a list
with open('output/AI/alleleCountsplits.txt') as f:
    splitList = f.read().splitlines()

numSplits = math.ceil(len(splitList)/2000)
splitListSubset = np.array_split(splitList, numSplits)

splitListSubset_dict = {}
for index, a in enumerate(splitListSubset):
    splitListSubset_dict[str(index)] = a

rule all:
    input:
        [expand('output/AI/alleleCountSplits/{splitFile}_' + str(config['minHets']) + 'hets', splitFile = l) for l in splitList],
        [expand('output/AI/alleleCounts_' + str(config['minHets']) + 'hets{splitGroup}.csv', splitGroup = key) for key in splitListSubset_dict],
        'output/AI/alleleCounts_' + str(config['minHets']) + 'hets.csv'

rule filterVariantHets:
    input:
        lambda wildcards: ['output/AI/alleleCountSplits/{splitFile}'.format(splitFile=wildcards.splitFile)],
        numVariantHets = config['numVariantHets']
    output:
        'output/AI/alleleCountSplits/{splitFile}_' + str(config['minHets']) + 'hets'
    params:
        minHets = config['minHets']
    log:
        out = 'output/AI/logs/{splitFile}_filterVariantHets.out'
    shell:
        """
        module load python/3.9.6
        python3 scripts/filterVariantHets.py {input} {params.minHets} 1> {log.out}
        """

rule concatVariantHets_part1:
    input:
        lambda wildcards: expand('output/AI/alleleCountSplits/{splitSubset}_' + str(config['minHets']) + 'hets', splitSubset = splitListSubset_dict[wildcards.splitGroup])
    output:
        'output/AI/alleleCounts_' + str(config['minHets']) + 'hets{splitGroup}.csv'
    params:
        minHets = config['minHets'],
        splitGroup = lambda wildcards: wildcards.splitGroup
    log:
        out = 'output/AI/logs/{splitGroup}_concatVariantHets.out'
    shell:
        """
        module load python/3.9.6
        python3 scripts/concatVariantHets_part1.py {params.minHets} {params.splitGroup} {input} 1> {log.out}
        """

rule concatVariantHets_part2:
    input:
        [expand("output/AI/alleleCounts_" + str(config['minHets']) + 'hets{splitGroup}.csv', splitGroup = key) for key in splitListSubset_dict]
    output:
        'output/AI/alleleCounts_' + str(config['minHets']) + 'hets.csv'
    params:
        minHets = config['minHets']
    log:
        out = 'output/AI/logs/concatVariantHets_part2.out'
    shell:
        """
        module load python/3.9.6
        python3 scripts/concatVariantHets_part2.py {params.minHets} {input} 1> {log.out}
        """



