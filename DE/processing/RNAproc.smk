#!/usr/bin/env python3

import pandas as pd
import os, shutil
import glob
from datetime import date

## Load config file
configfile: "config/config.yaml"

## Read in samplesheet
samples = pd.read_csv(config["samplesheet"],sep = ",")

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

include: "../../rules/RNAprocessing.smk"

rule all:
    input:
        'output/qc/multiqc_report.html',
        'output/' + str(date.today()) + '_gse.rda'

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

rule multiqc:
    input:
        [expand('output/qc/{group}_{R}_fastqc.zip', group = key, R = ['R1', 'R2']) for key in read1],
        [expand('output/{group}/trim/{group}_{R}.fastq.gz_trimming_report.txt', group = key, R =['R1', 'R2']) for key in read1],
        [expand('output/quant/{group}/quant.sf', group = key) for key in read1]
    output:
        'output/qc/multiqc_report.html'
    params:
        version = config['multiqcVersion']

    log:
        out = 'output/logs/multiqc.out',
        err = 'output/logs/multiqc.err'
    shell:
        """
        module load multiqc/{params.version}
        multiqc -f -o output/qc 1> {log.out} 2> {log.err}
        """

rule tximport:
    input:
        [expand('output/quant/{group}/quant.sf', group = key) for key in read1]
    output:
        'output/' + str(date.today()) + '_gse.rda'
    params:
        version = config['rVersion'],
        samplesheet = samples
    log:
        out = 'output/logs/tximport.out',
        err = 'output/logs/tximport.err'
    shell:
        """
        module load r/{params.version}
        Rscript tximport.R {params.samplesheet} 1> {log.out} 2> {log.err}
        """