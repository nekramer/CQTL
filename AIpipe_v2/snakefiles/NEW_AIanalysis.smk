#!/usr/bin/env python3

import pandas as pd
import os, re

## Load config file
configfile: "config/config_AIanalysis.yaml"

## Read in RNA-seq samplesheet
samples = pd.read_csv(config["samplesheet"],sep=",")

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

## Extract grouped donors
groupedDonors = samples.groupby('Donor')['Sample'].apply(list).to_dict()

## Get vcf file path and prefix of VCFproc processed vcf from config file
vcf = config['vcf']
vcf_file = os.path.basename(vcf)
vcf_prefix = vcf_file[:re.search("_nodups_biallelic.vcf.gz", vcf_file).span()[0]]


rule all:
    input:
        'output/AI/variants.csv',
        [expand("output/{group}/alleleCounts/{group}_alleleCounts_joined.csv", group = key) for key in read1],
        'output/AI/alleleCounts.pkl',
        [expand("output/AI/{donor}_alleleCounts_joined.csv", donor = key) for key in groupedDonors]
        #'output/AI/genohets.pkl'


rule getVariants:
    input: 
        [expand("output/{group}/alleleCounts/{group}_alleleCounts.csv", group = key) for key in read1]
    output:
        'output/AI/variants.csv'
    log:
        out = 'output/AI/logs/getVariants.out'
    shell:
        """
        module load r/4.2.1
        mkdir -p output/AI/logs
        Rscript scripts/getVariants.r {input} 1> {log.out}
        """

rule mergeSampleVariants:
    input:
        lambda wildcards: ['output/{group}/alleleCounts/{group}_alleleCounts.csv'.format(group=wildcards.group)],
        variants = rules.getVariants.output
    output:
        'output/{group}/alleleCounts/{group}_alleleCounts_joined.csv'
    params:
        minTotalAlleleCounts = config['minTotalAlleleCounts'],
        minAlleleCounts = config['minAlleleCounts']
    log:
        out = 'output/AI/logs/{group}_mergeSampleVariants.out'
    shell:
        """
        module load python/3.9.6
        mkdir -p output/AI/logs
        python3 scripts/mergeSampleVariants.py {input} {params.minTotalAlleleCounts} {params.minAlleleCounts} 1> {log.out}
        """

rule concatDonorConditions:
    input:
        lambda wildcards: expand('output/{group}/alleleCounts/{group}_alleleCounts_joined.csv', group = groupedDonors[wildcards.donor])
    output:
        'output/AI/{donor}_alleleCounts_joined.csv'
    params:
        donor = lambda wildcards: wildcards.donor
    log:
        out = 'output/AI/logs/{donor}_concatDonorConditions.out'
    shell:
        """
        module load python/3.9.6
        mkdir -p output/AI/logs
        python3 scripts/concatDonorConditions.py {params.donor} {input} 1> {log.out}
        """

rule getGenoHets:
    input:
        vcf = 'output/vcf/' + vcf_prefix + '_nodups_biallelic.vcf.gz',
        index = 'output/vcf/' + vcf_prefix + '_nodups_biallelic.vcf.gz.tbi'
    output:
        'output/AI/genohets.csv'
    params:
        donorConversions = 'donors.txt'
    log:
        out = "output/AI/logs/genoHets.out"
    shell:
        """
        python3 scripts/genoHets.py {input.vcf} {params.donorConversions} 1> {log.out}
        """

        

# rule checkDonorVariants:
#     input: 
#         rules.getGenoHets.output,
#         lambda wildcards: ['output/{group}/alleleCounts/{group}_alleleCounts_joined.csv'.format(group=wildcards.group)]
#     output:

#     log:

#     shell:


     
rule concatAlleleCounts:
    input: 
        [expand("output/{group}/alleleCounts/{group}_alleleCounts_joined.csv", group = key) for key in read1]
    output:
        'output/AI/alleleCounts.pkl'
    log:
        out = 'output/AI/logs/concatAlleleCounts.out'
    shell:
        """
        module load python/3.9.6
        mkdir -p output/AI/logs
        python3 scripts/concatAlleleCounts.py {input} 1> {log.out}
        """


