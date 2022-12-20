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

rule_all_inputs = [[expand('output/normquant/{condition}_CPMadjTMM_invNorm.bed.gz', condition = ['ALL', 'CTL', 'FNF'])],
                    [expand('output/normquant/{condition}_CPMadjTMM_invNorm.bed.gz.tbi', condition = ['ALL', 'CTL', 'FNF'])],
                    [expand('output/covar/{condition}_PEERfactors_k{Nk}.txt', condition = ['CTL', 'FNF'], Nk = Nk)],
                    [expand('output/covar/{condition}_PEERfactors_k{Nk}_variance.txt', condition = ['CTL', 'FNF'], Nk = Nk)]]

        
include: 'genoCovariate.smk'
batches = ['RNAKitBatch', 'RNASequencingBatch', 'genoBatch', 'DNAKitBatch']


if config['iteratePEER'] == 'TRUE':
    peerCov = 'output/covar/{condition}_PEER_k{Nk}_genoPC'
    peerQTL_perm = 'output/qtl/{condition}_PEER_k{Nk}_genoPC'
    peerQTL_nominal = 'output/qtl/{condition}_PEER_k{Nk}_genoPC'
    nomThreshold = 'output/qtl/{condition}_PEER_k{Nk}_genoPC'
    nomFinal = 'output/qtl/{condition}_PEER_k{Nk}_genoPC'
    peerMultipleTesting_perm = 'output/qtl/{condition}_PEER_k{Nk}_genoPC'

    PEER_eQTL_out = 'output/logs/{condition}_PEERk{Nk}_genoPC'
    PEER_eQTL_err = 'output/logs/{condition}_PEERk{Nk}_genoPC'
    PEER_nominal_eQTL_out = 'output/logs/{condition}_PEERk{Nk}_genoPC'
    PEER_nominal_eQTL_err = 'output/logs/{condition}_PEERk{Nk}_genoPC'
    get_nomThreshold_out = 'output/logs/{condition}_PEERk{Nk}_genoPC'
    get_nomThreshold_err = 'output/logs/{condition}_PEERk{Nk}_genoPC'
    nomFilter_out = 'output/logs/{condition}_PEERk{Nk}_genoPC'
    nomFilter_err = 'output/logs/{condition}_PEERk{Nk}_genoPC'
    PEER_multipleTesting_perm_out = 'output/logs/{condition}_PEERk{Nk}_genoPC'
    PEER_multipleTesting_perm_err = 'output/logs/{condition}_PEERk{Nk}_genoPC'

else:
    peerCov = 'output/covar/{condition}_PEER_kneedle_genoPC'
    peerQTL_perm = 'output/qtl/{condition}_PEER_kneedle_genoPC'
    peerQTL_nominal = 'output/qtl/{condition}_PEER_kneedle_genoPC'
    nomThreshold = 'output/qtl/{condition}_PEER_kneedle_genoPC'
    nomFinal = 'output/qtl/{condition}_PEER_kneedle_genoPC'
    peerMultipleTesting_perm = 'output/qtl/{condition}_PEER_kneedle_genoPC'

    PEER_eQTL_out = 'output/logs/{condition}_PEERkneedle_eQTL.out'
    PEER_eQTL_err = 'output/logs/{condition}_PEERkneedle_eQTL.err'
    PEER_nominal_eQTL_out = 'output/logs/{condition}_PEERkneedle_genoPC'
    PEER_nominal_eQTL_err = 'output/logs/{condition}_PEERkneedle_genoPC'
    get_nomThreshold_out = 'output/logs/{condition}_PEERkneedle_genoPC'
    get_nomThreshold_err = 'output/logs/{condition}_PEERkneedle_genoPC'
    nomFilter_out = 'output/logs/{condition}_PEERkneedle_genoPC'
    nomFilter_err = 'output/logs/{condition}_PEERkneedle_genoPC'
    PEER_multipleTesting_perm_out = 'output/logs/{condition}_PEERkneedle_genoPC'
    PEER_multipleTesting_perm_err = 'output/logs/{condition}_PEERkneedle_genoPC'


for b in batches:
    b_include = config[b]
    if b_include == "TRUE":

        peerCov += '_{}'.format(b)
        peerQTL_perm += '_{}'.format(b)
        peerQTL_nominal += '_{}'.format(b)
        nomThreshold += '_{}'.format(b)
        nomFinal += '_{}'.format(b)
        peerMultipleTesting_perm += '_{}'.format(b)

        PEER_eQTL_out += '_{}'.format(b)
        PEER_eQTL_err += '_{}'.format(b)
        PEER_nominal_eQTL_out += '_{}'.format(b)
        PEER_nominal_eQTL_err += '_{}'.format(b)
        get_nomThreshold_out += '_{}'.format(b)
        get_nomThreshold_err += '_{}'.format(b)
        nomFilter_out += '_{}'.format(b)
        nomFilter_err += '_{}'.format(b)
        PEER_multipleTesting_perm_out += '_{}'.format(b)
        PEER_multipleTesting_perm_err += '_{}'.format(b)

peerCov += ".txt"
peerQTL_perm += "_perm1Mb.txt"
peerQTL_nominal += "_nom1Mb.txt"
nomThreshold += "_nom1Mb_thresholds.csv"
nomFinal += "_nom1Mb_final.txt"
peerMultipleTesting_perm += "_perm1Mb_FDR.txt"

PEER_eQTL_out += "_eQTL.out"
PEER_eQTL_err += "_eQTL.err"
PEER_nominal_eQTL_out += "_nominal_eQTL.out"
PEER_nominal_eQTL_err += "_nominal_eQTL.err"
get_nomThreshold_out += "_nomThreshold.out"
get_nomThreshold_err += "_nomThreshold.err"
nomFilter_out += '_nomFilter.out'
nomFilter_err += '_nomFilter.err'
PEER_multipleTesting_perm_out += "_eQTL_multipleTesting_perm.out"
PEER_multipleTesting_perm_err += "_eQTL_multipleTesting_perm.err"


if config['iteratePEER'] == 'TRUE':
    rule_all_inputs.extend([[expand(peerCov, condition = ['CTL', 'FNF'], Nk = n) for n in range(5, Nk + 1, config['iterateBy'])]])
    rule_all_inputs.extend([[expand(peerQTL_perm, condition = ['CTL', 'FNF'], Nk = n) for n in range(5, Nk + 1, config['iterateBy'])]])
    rule_all_inputs.extend([[expand(peerQTL_nominal, condition = ['CTL', 'FNF'], Nk = n) for n in range(5, Nk + 1, config['iterateBy'])]])
    rule_all_inputs.extend([[expand(nomThreshold, condition = ['CTL', 'FNF'], Nk = n) for n in range(5, Nk + 1, config['iterateBy'])]])
    rule_all_inputs.extend([[expand(nomFinal, condition = ['CTL', 'FNF'], Nk = n) for n in range(5, Nk + 1, config['iterateBy'])]])
    rule_all_inputs.extend([[expand(peerMultipleTesting_perm, condition = ['CTL', 'FNF'], Nk = n) for n in range(5, Nk + 1, config['iterateBy'])]])
else:
    rule_all_inputs.extend([[expand(peerCov, condition = ['CTL', 'FNF'])]])
    rule_all_inputs.extend([[expand(peerQTL_perm, condition = ['CTL', 'FNF'])]])
    rule_all_inputs.extend([[expand(peerQTL_nominal, condition = ['CTL', 'FNF'])]])
    rule_all_inputs.extend([[expand(nomThreshold, condition = ['CTL', 'FNF'])]])
    rule_all_inputs.extend([[expand(nomFinal, condition = ['CTL', 'FNF'])]])
    rule_all_inputs.extend([[expand(peerMultipleTesting_perm, condition = ['CTL', 'FNF'])]])
    rule_all_inputs.extend([[expand('output/covar/{condition}_PEERkneedle.txt', condition = ['CTL', 'FNF'])]])

## Define rules
rule all:
    input:
        rule_all_inputs

include: "eQTL.smk"

# Permutation pass
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
        out = PEER_eQTL_out,
        err = PEER_eQTL_err
    shell:
        """
        module load qtltools/{params.version}
        QTLtools cis --vcf {input.vcf} --bed {input.bed} --cov {input.cov} --permute 1000 --window 1000000 --out {output} 1> {log.out} 2> {log.err}
        """

# Nominal pass    
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
        out = PEER_nominal_eQTL_out,
        err = PEER_nominal_eQTL_err
    shell:
        """
        module load qtltools/{params.version}
        QTLtools cis --vcf {input.vcf} --bed {input.bed} --cov {input.cov} --nominal 1 --window 1000000 --out {output} 1> {log.out} 2> {log.err}
        """ 

# Determining other significant variants using threshold determined by lead variants in permutation pass
rule get_nomThreshold:
    input:
        permData = rules.PEER_eQTL.output
    output:
        nomThreshold
    params:
        version = config['Rversion'],
        FDRthreshold = config['FDRthreshold']
    log:
        out = get_nomThreshold_out,
        err = get_nomThreshold_err
    shell:
        """
        module load r/{params.version}
        Rscript scripts/get_nomThreshold.R {input.permData} {params.FDRthreshold} {output} 1> {log.out} 2> {log.err}
        """

# Filtering nominal variants based on nominal threshold
rule nominal_Filter:
    input:
        nomData = rules.PEER_nominal_eQTL.output,
        nomThreshold = rules.get_nomThreshold.output
    output:
        nomFinal
    params:
        version = config['pythonVersion']
    log:
        out = nomFilter_out,
        err = nomFilter_err
    shell:
        """
        module load python/{params.version}
        python3 scripts/correct_nomQTLs.py {input.nomData} {input.nomThreshold} {output} 1> {log.out} 2> {log.err}
        """

# Add multiple testing correction to permutation pass lead variants/eGenes 
rule PEER_multipleTesting_perm:
    input:
        qtlResult = peerQTL_perm,
        geneInfo = rules.indexQuant.output.bed
    output:
        peerMultipleTesting_perm
    params:
        version = config['Rversion']
    log:
        out = PEER_multipleTesting_perm_out,
        err = PEER_multipleTesting_perm_err,
    shell:
        """
        module load r/{params.version}
        Rscript scripts/correctQTLs.R {input.qtlResult} {input.geneInfo} {output} 1> {log.out} 2> {log.err}
        """