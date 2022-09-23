#!/usr/bin/env python3
import pandas as pd
import os, shutil
import re
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
vcf_file = os.path.basename(vcf)
vcf_prefix = vcf_file[:re.search("_ALL_qc.vcf.gz", vcf_file).span()[0]]

## Define rules
rule all:
    input:
        [expand("output/mbv/{group}.bamstat.txt", group = key) for key in read1],
        [expand("output/quant/{group}/quant.sf", group = key) for key in read1]


include: "../../rules/VCFprocessing.smk"

include: "../../rules/RNAprocessing.smk"

rule mbv:
    input:
        bam = rules.align.output,
        index = rules.index.output,
        vcf = rules.zipVCF2.output
    output:
        'output/mbv/{group}.bamstat.txt'
    params:
        version = config['QTLToolsVersion']
    shell:
        """
        module load qtltools/{params.version}
        QTLtools mbv --bam {input.bam} --vcf {input.vcf} --out output/mbv/all.bamstat.txt
        """

# rule mbvParse:

rule quant:
    input:
        trim1 = rules.trim.output.trim1,
        trim2 = rules.trim.output.trim2
    output:
        "output/quant/{group}/quant.sf"
    params:
        version = config['salmonVersion'],
        index = config['salmon']
    log:
        out = 'output/logs/quant_{group}.out',
        err = 'output/logs/quant_{group}.err'

    shell:
        """
        module load salmon/{params.version}
        salmon quant --writeUnmappedNames -l A -1 {input.trim1} -2 {input.trim2} -i {params.index} -o output/quant/{wildcards.group} --threads 1 1> {log.out} 2> {log.err}
        """

