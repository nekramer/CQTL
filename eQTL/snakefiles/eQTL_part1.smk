#!/usr/bin/env python3

import pandas as pd
import os, shutil
import glob

## Load config file
configfile: "config/config.yaml"

## Read in samplesheet
samples = pd.read_csv(config["samplesheet"], sep = ",")

## Convert samplesheet columns to strings
samples = samples.astype(str)

## Concatenate Sequencing_Directory to Read1 and Read2 for full read paths
samples['Read1'] = samples[['Sequencing_Directory', 'Read1']].apply(lambda row: os.path.join(*row), axis=1)
samples['Read2'] = samples[['Sequencing_Directory', 'Read2']].apply(lambda row: os.path.join(*row), axis=1)

## Group Seq_Reps
samples['id'] = samples[['Proj', 'Donor']].agg('_'.join, axis=1) + '_R_' + samples[['Condition', 'Time', 'Tech_Rep']].agg('_'.join, axis=1)

## Extract grouped read1 and read2s
read1 = samples.groupby(['id'])['Read1'].apply(list).to_dict()
read2 = samples.groupby(['id'])['Read2'].apply(list).to_dict()

## Get vcf file path of post-imputed, qc'd gzipped vcf file
vcf = config["vcf"]

## Define rules
rule all:
    input:
        [expand("output/mbv/{group}.bamstat.txt", group = key) for key in read1]


include: "../../rules/catR1.rule"

include: "../../rules/catR2.rule"

include: "../../rules/qc.rule"

include: "../../rules/trim.rule"

include: "../../rules/align.rule"

rule mbv:
    input:
        bam = rules.align.output,
        vcf = vcf
    output:
        'output/mbv/{group}.bamstat.txt'
    params:
        version = config['QTLToolsVersion']
    shell:
        """
        module load qtltools/{params.version}
        QTLtools mbv --bam {input.bam} --vcf {input.vcf} --out output/mbv/{wildcards.group}.bamstat.txt
        """




