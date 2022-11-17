#!/usr/bin/env python3
import pandas as pd
import os, shutil
import re
import glob
import numpy as np
from kneed import KneeLocator

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

rule_all_inputs = [[expand("output/covar/{condition}_PC.csv", condition = ['CTL', 'FNF', 'ALL'])],
        [expand('output/pca/{condition}_PC.csv', condition = ['CTL', 'FNF'])],
        [expand('output/covar/{condition}_PCkneedle.txt', condition = ['CTL', 'FNF'])],
        [expand('output/covar/{condition}_PEERfactors_k{Nk}.txt', condition = ['CTL', 'FNF'], Nk = n) for n in range(1, 16, 1)],
        [expand('output/plots/{condition}_PEER{Nk}_correlation.png', condition = ['CTL', 'FNF'], Nk = n) for n in range(1, 16, 1)],
        [expand('output/covar/{condition}_validPEER.txt', condition = ['CTL', 'FNF'])],
        [expand('output/plots/{condition}_pca.pdf', condition = ['CTL', 'FNF'])],
        [expand('output/plots/{condition}_screeplot.pdf', condition = ['CTL', 'FNF'])],
        [expand('output/plots/{condition}_PC_correlation.png', condition = ['CTL', 'FNF'])]]


        #[expand("output/mbv/{group}.bamstat.txt", group = key) for key in read1]

if config['genoCovar'] == 'yes':
    include: 'genoCovariate.smk'
    if config['batchCovar'] == 'TRUE':
        rule_all_inputs.extend([[expand('output/covar/{condition}_PC_genoPC_batch_covar.txt', condition = ['CTL', 'FNF'])]])
        pcCov = ['output/covar/{condition}_PC_genoPC_batch_covar.txt']

        rule_all_inputs.extend([[expand('output/covar/{condition}_PEER_k{Nk}_genoPC_batch_covar.txt', condition = ['CTL', 'FNF'], Nk = n) for n in range(1, 16, 1)]])
        peerCov = ['output/covar/{condition}_PEER_k{Nk}_genoPC_batch_covar.txt']

        
        rule_all_inputs.extend([[expand('output/qtl/{condition}_PC_genoPC_batch_perm1Mb.txt', condition = ['CTL', 'FNF'])]])
        pcQTL_perm = ['output/qtl/{condition}_PC_genoPC_batch_perm1Mb.txt']
        
        rule_all_inputs.extend([[expand('output/qtl/{condition}_PEER_k{Nk}_genoPC_batch_perm1Mb.txt', condition = ['CTL', 'FNF'], Nk = n) for n in range(1, 16, 1)]])
        peerQTL_perm = ['output/qtl/{condition}_PEER_k{Nk}_genoPC_batch_perm1Mb.txt']
        

        rule_all_inputs.extend([[expand('output/qtl/{condition}_PC_genoPC_batch_nom1Mb.txt', condition = ['CTL', 'FNF'])]])
        pcQTL_nominal = ['output/qtl/{condition}_PC_genoPC_batch_nom1Mb.txt']

        rule_all_inputs.extend([[expand('output/qtl/{condition}_PEER_k{Nk}_genoPC_batch_nom1Mb.txt', condition = ['CTL', 'FNF'], Nk = n) for n in range(1, 16, 1)]])
        peerQTL_nominal = ['output/qtl/{condition}_PEER_k{Nk}_genoPC_batch_nom1Mb.txt']
        

        rule_all_inputs.extend([[expand('output/qtl/{condition}_PEER_k{Nk}_genoPC_batch_nom1Mb_final.txt', condition = ['CTL', 'FNF'], Nk = n) for n in range(1, 16, 1)]])
        PEER_postProcessing = ['output/qtl/{condition}_PEER_k{Nk}_genoPC_batch_nom1Mb_final.txt']

        rule_all_inputs.extend([[expand('output/qtl/{condition}_PC_genoPC_batch_nom1Mb_final.txt', condition = ['CTL', 'FNF'])]])
        PC_postProcessing = ['output/qtl/{condition}_PC_genoPC_batch_nom1Mb_final.txt']

        rule_all_inputs.extend([[expand('output/qtl/{condition}_PC_genoPC_batch_perm1Mb_FDR.txt', condition = ['CTL', 'FNF'])]])
        pcMultipleTesting_perm = ['output/qtl/{condition}_PC_genoPC_batch_perm1Mb_FDR.txt']

        rule_all_inputs.extend([[expand('output/qtl/{condition}_PEER_k{Nk}_genoPC_batch_perm1Mb_FDR.txt', condition = ['CTL', 'FNF'], Nk = n) for n in range(1, 16, 1)]])
        peerMultipleTesting_perm = ['output/qtl/{condition}_PEER_k{Nk}_genoPC_batch_perm1Mb_FDR.txt']


    else:
        rule_all_inputs.extend([[expand('output/covar/{condition}_PC_genoPC_covar.txt', condition = ['CTL', 'FNF'])]])
        pcCov = ['output/covar/{condition}_PC_genoPC_covar.txt']

        rule_all_inputs.extend([[expand('output/covar/{condition}_PEER_k{Nk}_genoPC_covar.txt', condition = ['CTL', 'FNF'], Nk = n) for n in range(1, 16, 1)]])
        peerCov = ['output/covar/{condition}_PEER_k{Nk}_genoPC_covar.txt']
        
        rule_all_inputs.extend([[expand('output/qtl/{condition}_PC_genoPC_perm1Mb.txt', condition = ['CTL', 'FNF'])]])
        pcQTL_perm = ['output/qtl/{condition}_PC_genoPC_perm1Mb.txt']

        rule_all_inputs.extend([[expand('output/qtl/{condition}_PEER_k{Nk}_genoPC_perm1Mb.txt', condition = ['CTL', 'FNF'], Nk = n) for n in range(1, 16, 1)]])
        peerQTL_perm = ['output/qtl/{condition}_PEER_k{Nk}_genoPC_perm1Mb.txt']

        rule_all_inputs.extend([[expand('output/qtl/{condition}_PC_genoPC_nom1Mb.txt', condition = ['CTL', 'FNF'])]])
        pcQTL_nominal = ['output/qtl/{condition}_PC_genoPC_nom1Mb.txt']

        rule_all_inputs.extend([[expand('output/qtl/{condition}_PEER_k{Nk}_genoPC_nom1Mb.txt', condition = ['CTL', 'FNF'], Nk = n) for n in range(1, 16, 1)]])
        peerQTL_nominal = ['output/qtl/{condition}_PEER_k{Nk}_genoPC_nom1Mb.txt']

        rule_all_inputs.extend([[expand('output/qtl/{condition}_PC_genoPC_nom1Mb_final.txt', condition = ['CTL', 'FNF'])]])
        PC_postProcessing = ['output/qtl/{condition}_PC_genoPC_nom1Mb_final.txt']

        rule_all_inputs.extend([[expand('output/qtl/{condition}_PEER_k{Nk}_genoPC_nom1Mb_final.txt', condition = ['CTL', 'FNF'], Nk = n) for n in range(1, 16, 1)]])
        PEER_postProcessing = ['output/qtl/{condition}_PEER_k{Nk}_genoPC_nom1Mb_final.txt']

        rule_all_inputs.extend([[expand('output/qtl/{condition}_PC_genoPC_perm1Mb_FDR.txt', condition = ['CTL', 'FNF'])]])
        pcMultipleTesting_perm = ['output/qtl/{condition}_PC_genoPC_perm1Mb_FDR.txt']

        rule_all_inputs.extend([[expand('output/qtl/{condition}_PEER_k{Nk}_genoPC_perm1Mb_FDR.txt', condition = ['CTL', 'FNF'], Nk = n) for n in range(1, 16, 1)]])
        peerMultipleTesting_perm = ['output/qtl/{condition}_PEER_k{Nk}_genoPC_perm1Mb_FDR.txt']

else:
    include: 'no_genoCovariate.smk'
    if config['batchCovar'] == 'TRUE':
        rule_all_inputs.extend([[expand('output/covar/{condition}_PC_batch_covar.txt', condition = ['CTL', 'FNF'])]])
        pcCov = ['output/covar/{condition}_PC_batch_covar.txt']

        rule_all_inputs.extend([[expand('output/covar/{condition}_PEER_k{Nk}_batch_covar.txt', condition = ['CTL', 'FNF'], Nk = n) for n in range(1, 16, 1)]])
        peerCov = ['output/covar/{condition}_PEER_k{Nk}_batch_covar.txt']

        rule_all_inputs.extend([[expand('output/qtl/{condition}_PC_batch_perm1Mb.txt', condition = ['CTL', 'FNF'])]])
        pcQTL_perm = ['output/qtl/{condition}_PC_batch_perm1Mb.txt']

        rule_all_inputs.extend([[expand('output/qtl/{condition}_PEER_k{Nk}_batch_perm1Mb.txt', condition = ['CTL', 'FNF'], Nk = n) for n in range(1, 16, 1)]])
        peerQTL_perm = ['output/qtl/{condition}_PEER_k{Nk}_batch_perm1Mb.txt']
        
        rule_all_inputs.extend([[expand('output/qtl/{condition}_PC_batch_nom1Mb.txt', condition = ['CTL', 'FNF'])]])
        pcQTL_nominal = ['output/qtl/{condition}_PC_batch_nom1Mb.txt']

        rule_all_inputs.extend([[expand('output/qtl/{condition}_PEER_k{Nk}_batch_nom1Mb.txt', condition = ['CTL', 'FNF'], Nk = n) for n in range(1, 16, 1)]])
        peerQTL_nominal = ['output/qtl/{condition}_PEER_k{Nk}_batch_nom1Mb.txt']

        rule_all_inputs.extend([[expand('output/qtl/{condition}_PC_batch_nom1Mb_final.txt', condition = ['CTL', 'FNF'])]])
        PC_postProcessing = ['output/qtl/{condition}_PC_batch_nom1Mb_final.txt']

        rule_all_inputs.extend([[expand('output/qtl/{condition}_PEER_k{Nk}_batch_nom1Mb_final.txt', condition = ['CTL', 'FNF'], Nk = n) for n in range(1, 16, 1)]])
        PEER_postProcessing = ['output/qtl/{condition}_PEER_k{Nk}_batch_nom1Mb_final.txt']

        rule_all_inputs.extend([[expand('output/qtl/{condition}_PC_batch_perm1Mb_FDR.txt', condition = ['CTL', 'FNF'])]])
        pcMultipleTesting_perm = ['output/qtl/{condition}_PC_batch_perm1Mb_FDR.txt']
        
        rule_all_inputs.extend([[expand('output/qtl/{condition}_PEER_k{Nk}_batch_perm1Mb_FDR.txt', condition = ['CTL', 'FNF'], Nk = n) for n in range(1, 16, 1)]])
        peerMultipleTesting_perm = ['output/qtl/{condition}_PEER_k{Nk}_batch_perm1Mb_FDR.txt']

    else:
        rule_all_inputs.extend([[expand('output/covar/{condition}_PC_covar.txt', condition = ['CTL', 'FNF'])]])
        pcCov = ['output/covar/{condition}_PC_covar.txt']

        rule_all_inputs.extend([[expand('output/covar/{condition}_PEER_k{Nk}_covar.txt', condition = ['CTL', 'FNF'], Nk = n) for n in range(1, 16, 1)]])
        peerCov = ['output/covar/{condition}_PEER_k{Nk}_covar.txt']
        
        rule_all_inputs.extend([[expand('output/qtl/{condition}_PC_perm1Mb.txt', condition = ['CTL', 'FNF'])]])
        pcQTL_perm = ['output/qtl/{condition}_PC_perm1Mb.txt']

        rule_all_inputs.extend([[expand('output/qtl/{condition}_PEER_k{Nk}_perm1Mb.txt', condition = ['CTL', 'FNF'], Nk = n) for n in range(1, 16, 1)]])
        peerQTL_perm = ['output/qtl/{condition}_PEER_k{Nk}_perm1Mb.txt']
        
        rule_all_inputs.extend([[expand('output/qtl/{condition}_PC_nom1Mb.txt', condition = ['CTL', 'FNF'])]])
        pcQTL_nominal = ['output/qtl/{condition}_PC_nom1Mb.txt']

        rule_all_inputs.extend([[expand('output/qtl/{condition}_PEER_k{Nk}_nom1Mb.txt', condition = ['CTL', 'FNF'], Nk = n) for n in range(1, 16, 1)]])
        peerQTL_nominal = ['output/qtl/{condition}_PEER_k{Nk}_nom1Mb.txt']

        rule_all_inputs.extend([[expand('output/qtl/{condition}_PC_nom1Mb_final.txt', condition = ['CTL', 'FNF'])]])
        PC_postProcessing = ['output/qtl/{condition}_PC_nom1Mb_final.txt']

        rule_all_inputs.extend([[expand('output/qtl/{condition}_PEER_k{Nk}_nom1Mb_final.txt', condition = ['CTL', 'FNF'], Nk = n) for n in range(1, 16, 1)]])
        PEER_postProcessing = ['output/qtl/{condition}_PEER_k{Nk}_nom1Mb_final.txt']

        rule_all_inputs.extend([[expand('output/qtl/{condition}_PC_perm1Mb_FDR.txt', condition = ['CTL', 'FNF'])]])
        pcMultipleTesting_perm = ['output/qtl/{condition}_PC_perm1Mb_FDR.txt']
        
        rule_all_inputs.extend([[expand('output/qtl/{condition}_PEER_k{Nk}_perm1Mb_FDR.txt', condition = ['CTL', 'FNF'], Nk = n) for n in range(1, 16, 1)]])
        peerMultipleTesting_perm = ['output/qtl/{condition}_PEER_k{Nk}_perm1Mb_FDR.txt']

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

# rule mbvParse:

rule PC_eQTL:
    input:
        vcf = rules.renameVCFdonors.output.vcf,
        vcfIndex = rules.renameVCFdonors.output.index,
        bed = rules.indexQuant.output.bed,
        bedIndex = rules.indexQuant.output.index,
        cov = pcCov
    output:
        pcQTL_perm
    params:
        version = config['QTLToolsVersion']
    log:
        out = 'output/logs/{condition}_PCeQTL.out',
        err = 'output/logs/{condition}_PCeQTL.err'
    shell:
        """
        module load qtltools/{params.version}
        QTLtools cis --vcf {input.vcf} --bed {input.bed} --cov {input.cov} --permute 1000 --window 1000000 --out {output} 1> {log.out} 2> {log.err}
        """

rule PEER_eQTL:
    input:
        vcf = rules.renameVCFdonors.output.vcf,
        vcfIndex = rules.renameVCFdonors.output.index,
        bed = rules.indexQuant.output.bed,
        bedIndex = rules.indexQuant.output.index,
        cov = peerCov,
        check = lambda wildcards: 'output/covar/{condition}_validPEER.txt'.format(condition=wildcards.condition)
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
        if grep -Fx {wildcards.Nk} {input.check}
        then
            QTLtools cis --vcf {input.vcf} --bed {input.bed} --cov {input.cov} --permute 1000 --window 1000000 --out {output} 1> {log.out} 2> {log.err}
        else
            touch {output}
        fi
        """
    
rule PC_nominal_eQTL:
    input:
        vcf = rules.renameVCFdonors.output.vcf,
        vcfIndex = rules.renameVCFdonors.output.index,
        bed = rules.indexQuant.output.bed,
        bedIndex = rules.indexQuant.output.index,
        cov = pcCov
    output:
        pcQTL_nominal
    params:
        version = config['QTLToolsVersion']
    log:
        out = 'output/logs/{condition}_PCnominal_eQTL.out',
        err = 'output/logs/{condition}_PCnominal_eQTL.err'
    shell:
        """
        module load qtltools/{params.version}
        QTLtools cis --vcf {input.vcf} --bed {input.bed} --cov {input.cov} --nominal 0.01 --window 1000000 --out {output} 1> {log.out} 2> {log.err}
        """


rule PEER_nominal_eQTL:
    input:
        vcf = rules.renameVCFdonors.output.vcf,
        vcfIndex = rules.renameVCFdonors.output.index,
        bed = rules.indexQuant.output.bed,
        bedIndex = rules.indexQuant.output.index,
        cov = peerCov,
        check = lambda wildcards: 'output/covar/{condition}_validPEER.txt'.format(condition=wildcards.condition)
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
        if grep -Fx {wildcards.Nk} {input.check}
        then
            QTLtools cis --vcf {input.vcf} --bed {input.bed} --cov {input.cov} --nominal 1 --window 1000000 --out {output} 1> {log.out} 2> {log.err}
        else
            touch {output}
        fi
        """ 

rule PC_qtlPostProcessing:
    input:
        nomData = rules.PC_nominal_eQTL.output,
        permData = rules.PC_eQTL.output,
        geneInfo = rules.indexQuant.output.bed
    output:
        PC_postProcessing
    params:
        version = config['Rversion'],
        FDRthreshold = config['FDRthreshold']
    log:
        out = 'output/logs/{condition}_PC_qtlPostProcessing.out',
        err = 'output/logs/{condition}_PC_qtlPostProcessing.err',
    shell:
        """
        module load r/{params.version}
        Rscript scripts/postProcessQTLtools.R {input.nomData} {input.permData} {params.FDRthreshold} {input.geneInfo} {output} 1> {log.out} 2> {log.err}
        """       

rule PEER_qtlPostProcessing:
    input:
        nomData = rules.PEER_nominal_eQTL.output,
        permData = rules.PEER_eQTL.output,
        geneInfo = rules.indexQuant.output.bed,
        check = lambda wildcards: 'output/covar/{condition}_validPEER.txt'.format(condition=wildcards.condition)
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
        if grep -Fx {wildcards.Nk} {input.check}
        then
            Rscript scripts/postProcessQTLtools.R {input.nomData} {input.permData} {params.FDRthreshold} {input.geneInfo} {output} 1> {log.out} 2> {log.err}
        else
            touch {output}
        fi
        """

rule PC_multipleTesting_perm:
    input:
        qtlResult = pcQTL_perm,
        geneInfo = rules.indexQuant.output.bed
    output:
        pcMultipleTesting_perm
    params:
        version = config['Rversion']
    log:
        out = 'output/logs/{condition}_PCeQTL_multipleTesting_perm.out',
        err = 'output/logs/{condition}_PCeQTL_multipleTesting_perm.err'
    shell:
        """
        module load r/{params.version}
        Rscript scripts/correctQTLs.R {input.qtlResult} {input.geneInfo} {output} 1> {log.out} 2> {log.err}
        """

rule PEER_multipleTesting_perm:
    input:
        qtlResult = peerQTL_perm,
        geneInfo = rules.indexQuant.output.bed,
        check = lambda wildcards: 'output/covar/{condition}_validPEER.txt'.format(condition=wildcards.condition)
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
        if grep -Fx {wildcards.Nk} {input.check}
        then
            Rscript scripts/correctQTLs.R {input.qtlResult} {input.geneInfo} {output} 1> {log.out} 2> {log.err}
        else
            touch {output}
        fi
        """