#!/usr/bin/env python3
import os
import numpy as np
import math
import shutil

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

onsuccess:
    # Remove splits
    shutil.rmtree('output/AI/alleleCountSplits')
    os.remove('output/AI/alleleCountsplits.txt')

    # Remove noncombined alleleCount het files
    for key in splitListSubset_dict:
        os.remove('output/AI/alleleCounts_' + str(config['minHets']) + 'hets{splitGroup}.csv'.format(splitGroup = key))

rule all:
    input:
        [expand('output/AI/alleleCountSplits/{splitFile}_' + str(config['minHets']) + 'hets', splitFile = l) for l in splitList],
        [expand('output/AI/alleleCounts_' + str(config['minHets']) + 'hets{splitGroup}.csv', splitGroup = key) for key in splitListSubset_dict],
        'output/AI/alleleCounts_' + str(config['minHets']) + 'hets.csv',
        'output/AI/alleleCountsMatrix.csv',
        'output/AI/weightsMatrix.csv',
        'output/AI/colData.csv'

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

rule pivotMatrix:
    input:
        rules.concatVariantHets_part2.output
    output:
        'output/AI/alleleCountsMatrix.csv',
        'output/AI/weightsMatrix.csv',
        'output/AI/colData.csv'
    log:
        out = 'output/AI/logs/pivotMatrix.out'
    shell:
        """
        module load python/3.9.6
        python3 scripts/pivotMatrix.py {input} 1> {log.out}
        """

