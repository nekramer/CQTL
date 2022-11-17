#!/usr/bin/env python3
import pandas as pd
import os, shutil
import re
import glob
import numpy as np
from kneed import KneeLocator
import os.path
from os import path
import subprocess

## Load config file
configfile: "config/config_QTLtools_eQTL.yaml"

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

## Number of PEER factors
Nk = config['PEERfactors']

# rule_all_inputs = [[expand("output/covar/{condition}_PC.csv", condition = ['CTL', 'FNF', 'ALL'])],
#         [expand('output/pca/{condition}_PC.csv', condition = ['CTL', 'FNF'])],
#         [expand('output/covar/{condition}_PCkneedle.txt', condition = ['CTL', 'FNF'])],
#         [expand('output/covar/{condition}_PEERfactors_k{Nk}.txt', condition = ['CTL', 'FNF'], Nk = n) for n in range(1, 16, 1)],
#         [expand('output/plots/{condition}_PEER{Nk}_correlation.png', condition = ['CTL', 'FNF'], Nk = n) for n in range(1, 16, 1)],
#         [expand('output/covar/{condition}_validPEER.txt', condition = ['CTL', 'FNF'])],
#         [expand('output/plots/{condition}_pca.pdf', condition = ['CTL', 'FNF'])],
#         [expand('output/plots/{condition}_screeplot.pdf', condition = ['CTL', 'FNF'])],
#         [expand('output/plots/{condition}_PC_correlation.png', condition = ['CTL', 'FNF'])]]

rule_all_inputs = [[expand('output/covar/{condition}_PEERfactors_k{Nk}.txt', condition = ['CTL', 'FNF'], Nk = Nk)],
                    [expand('output/covar/{condition}_PEERfactors_k{Nk}_variance.txt', condition = ['CTL', 'FNF'], Nk = Nk)],
                    [expand('output/covar/{condition}_PEERkneedle.txt', condition = ['CTL', 'FNF'])]]

        #[expand("output/mbv/{group}.bamstat.txt", group = key) for key in read1]

include: 'genoCovariate.smk'
batches = ['RNAKitBatch', 'RNASequencingBatch', 'genoBatch', 'DNAKitBatch']


if config['iteratePEER'] == 'TRUE':
    peerCov = 'output/covar/{condition}_PEER_k{Nk}_genoPC'
    peerQTL_perm = 'output/qtl/{condition}_PEER_k{Nk}_genoPC'
    peerQTL_nominal = 'output/qtl/{condition}_PEER_k{Nk}_genoPC'
    PEER_postProcessing = 'output/qtl/{condition}_PEER_k{Nk}_genoPC'
    peerMultipleTesting_perm = 'output/qtl/{condition}_PEER_k{Nk}_genoPC'
else:
    peerCov = 'output/covar/{condition}_PEER_kneedle_genoPC'
    peerQTL_perm = 'output/qtl/{condition}_PEER_kneedle_genoPC'
    peerQTL_nominal = 'output/qtl/{condition}_PEER_kneedle_genoPC'
    PEER_postProcessing = 'output/qtl/{condition}_PEER_kneedle_genoPC'
    peerMultipleTesting_perm = 'output/qtl/{condition}_PEER_kneedle_genoPC'


for b in batches:
    b_include = config[b]
    if b_include == "TRUE":

        peerCov += '_{}'.format(b)
        peerQTL_perm += '_{}'.format(b)
        peerQTL_nominal += '_{}'.format(b)
        PEER_postProcessing += '_{}'.format(b)
        peerMultipleTesting_perm += '_{}'.format(b)

peerCov += ".txt"
peerQTL_perm += "_perm1Mb.txt"
peerQTL_nominal += "_nom1Mb.txt"
PEER_postProcessing += "_nom1Mb_final.txt"
peerMultipleTesting_perm += "_perm1Mb_FDR.txt"


if config['iteratePEER'] == 'TRUE':
    rule_all_inputs.extend([[expand(peerCov, condition = ['CTL', 'FNF'], Nk = n) for n in range(5, Nk + 1, 5)]])
    rule_all_inputs.extend([[expand(peerQTL_perm, condition = ['CTL', 'FNF'], Nk = n) for n in range(5, Nk + 1, 5)]])
    rule_all_inputs.extend([[expand(peerQTL_nominal, condition = ['CTL', 'FNF'], Nk = n) for n in range(5, Nk + 1, 5)]])
    rule_all_inputs.extend([[expand(PEER_postProcessing, condition = ['CTL', 'FNF'], Nk = n) for n in range(5, Nk + 1, 5)]])
    rule_all_inputs.extend([[expand(peerMultipleTesting_perm, condition = ['CTL', 'FNF'], Nk = n) for n in range(5, Nk + 1, 5)]])
else:
    rule_all_inputs.extend([[expand(peerCov, condition = ['CTL', 'FNF'])]])
    rule_all_inputs.extend([[expand(peerQTL_perm, condition = ['CTL', 'FNF'])]])
    rule_all_inputs.extend([[expand(peerQTL_nominal, condition = ['CTL', 'FNF'])]])
    rule_all_inputs.extend([[expand(PEER_postProcessing, condition = ['CTL', 'FNF'])]])
    rule_all_inputs.extend([[expand(peerMultipleTesting_perm, condition = ['CTL', 'FNF'])]])

## Define rules
rule all:
    input:
        rule_all_inputs

include: "eQTL.smk"

# rule mbv:
#     input:
#         bam = rules.align.output,
#         index = rules.index.output,
#         vcf = rules.renameVCFdonors.output.vcf
#     output:
#         'output/mbv/{group}.bamstat.txt'
#     params:
#         version = config['QTLToolsVersion']
#     shell:
#         """
#         module load qtltools/{params.version}
#         QTLtools mbv --bam {input.bam} --vcf {input.vcf} --out output/mbv/{wildcards.group}.bamstat.txt
#         """

# # rule mbvParse:

rule PEER_eQTL:
    input:
        vcf = rules.filterVCFvariants.output.vcf,
        vcfIndex = rules.filterVCFvariants.output.index,
        bed = rules.indexQuant.output.bed,
        bedIndex = rules.indexQuant.output.index,
        cov = peerCov
    output:
        peerQTL_perm
    params:
        version = config['QTLToolsVersion']
    log:
        out = 'output/logs/{condition}_PEERk{Nk}_eQTL.out',
        err = 'output/logs/{condition}_PEERk{Nk}_eQTL.err'
    shell:
        """
        module load qtltools/{params.version}
        QTLtools cis --vcf {input.vcf} --bed {input.bed} --cov {input.cov} --permute 1000 --window 1000000 --out {output} 1> {log.out} 2> {log.err}
        """
    
rule PEER_nominal_eQTL:
    input:
        vcf = rules.filterVCFvariants.output.vcf,
        vcfIndex = rules.filterVCFvariants.output.index,
        bed = rules.indexQuant.output.bed,
        bedIndex = rules.indexQuant.output.index,
        cov = peerCov
    output:
        peerQTL_nominal
    params:
        version = config['QTLToolsVersion']
    log:
        out = 'output/logs/{condition}_PEERk{Nk}_nominal_eQTL.out',
        err = 'output/logs/{condition}_PEERk{Nk}_nominal_eQTL.err'
    shell:
        """
        module load qtltools/{params.version}
        QTLtools cis --vcf {input.vcf} --bed {input.bed} --cov {input.cov} --nominal 1 --window 1000000 --out {output} 1> {log.out} 2> {log.err}
        """ 

rule PEER_qtlPostProcessing:
    input:
        nomData = rules.PEER_nominal_eQTL.output,
        permData = rules.PEER_eQTL.output,
        geneInfo = rules.indexQuant.output.bed
    output:
        PEER_postProcessing
    params:
        version = config['Rversion'],
        FDRthreshold = config['FDRthreshold']
    log:
        out = 'output/logs/{condition}_PEERk{Nk}_qtlPostProcessing.out',
        err = 'output/logs/{condition}_PEERk{Nk}_qtlPostProcessing.err',
    shell:
        """
        module load r/{params.version}
        Rscript scripts/postProcessQTLtools.R {input.nomData} {input.permData} {params.FDRthreshold} {input.geneInfo} {output} 1> {log.out} 2> {log.err}
        """

rule PEER_multipleTesting_perm:
    input:
        qtlResult = peerQTL_perm,
        geneInfo = rules.indexQuant.output.bed
    output:
        peerMultipleTesting_perm
    params:
        version = config['Rversion']
    log:
        out = 'output/logs/{condition}_PEERk{Nk}_eQTL_multipleTesting_perm.out',
        err = 'output/logs/{condition}_PEERk{Nk}_eQTL_multipleTesting_perm.err',
    shell:
        """
        module load r/{params.version}
        Rscript scripts/correctQTLs.R {input.qtlResult} {input.geneInfo} {output} 1> {log.out} 2> {log.err}
        """