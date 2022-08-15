#!/usr/bin/env python3

import pandas as pd
import os, re

## Load config file
configfile: "config/config_filterVariants.yaml"

## Get split file names in a list
with open('output/AI/alleleCountsplits.txt') as f:
    splitList = f.read().splitlines()

rule all:
    input:
        [expand('output/AI/alleleCountSplits/{splitFile}_' + str(config['minHets']) + 'hets', splitFile = l) for l in splitList],
        #'output/AI/alleleCounts_' + str(config['minHets']) + 'hets.csv'

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

# rule concatVariantHets:
#     input:
#          [expand("output/AI/alleleCountSplits/{splitFile}" + str(config['minHets']) + 'hets', splitFile = l) for l in splitList]
#     output:
#         'output/AI/alleleCounts_' + str(config['minHets']) + 'hets.csv'
#     params:
#         minHets = config['minHets']
#     log:
#         out = 'output/AI/logs/concatVariantHets.out'
#     shell:
#         """
#         module load python/3.9.6
#         python3 scripts/concatVariantHets.py {params.minHets} {input} 1> {log.out}
#         """


