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
        [expand("output/quant/{group}/quant.sf", group = key) for key in read1],
        [expand("output/normquant/{condition}_CPMadjTMM_invNorm.bed.gz", condition = ['CTL', 'FNF', 'ALL'])],
        [expand("output/normquant/{condition}_CPMadjTMM_invNorm.bed.gz.tbi", condition = ['CTL', 'FNF', 'ALL'])],
        [expand("output/covar/{condition}_PC.csv", condition = ['CTL', 'FNF', 'ALL'])],
        [expand('output/covar/{condition}_PCcovar.txt', condition = ['CTL', 'FNF', 'ALL'])],
        #[expand('output/qtl/PC_{condition}perm1Mb.txt', condition = ['CTL', 'FNF'])]


include: "../../rules/VCFprocessing.smk"

include: "../../rules/RNAprocessing.smk"

rule renameVCFdonors:
    input:
        #vcf = rules.zipVCF2.output
        vcf = rules.updateConfig.output.v
    output:
        vcf = 'output/vcf/' + vcf_prefix + '_newcontig_rename.vcf.gz',
        index = 'output/vcf/' + vcf_prefix + '_newcontig_rename.vcf.gz.tbi'
    params:
        donors = ",".join(samples['Donor'].unique().tolist()),
        samtoolsVersion = config['samtoolsVersion'],
        pythonVersion = config['pythonVersion'],
        prefix = vcf_prefix
    shell:
        """
        module load samtools/{params.samtoolsVersion}
        module load python/{params.pythonVersion}
        bcftools query -l {input.vcf} > donors.txt
        python3 scripts/renameVCFdonors.py donors.txt {params.donors}
        gunzip {input.vcf}
        bcftools reheader -s samples.txt -o output/vcf/{params.prefix}_newcontig_rename.vcf output/vcf/{params.prefix}_newcontig.vcf
        bgzip output/vcf/{params.prefix}_newcontig_rename.vcf && tabix -p vcf {output.vcf}
        """

rule mbv:
    input:
        bam = rules.align.output,
        index = rules.index.output,
        vcf = rules.renameVCFdonors.output.vcf
    output:
        'output/mbv/{group}.bamstat.txt'
    params:
        version = config['QTLToolsVersion']
    shell:
        """
        module load qtltools/{params.version}
        QTLtools mbv --bam {input.bam} --vcf {input.vcf} --out output/mbv/{wildcards.group}.bamstat.txt
        """

# # rule mbvParse:

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

rule quantNorm:
    input:
        [expand("output/quant/{group}/quant.sf", group = key) for key in read1]
    output:
        [expand("output/normquant/{condition}_CPMadjTMM_invNorm.bed", condition = ['CTL', 'FNF', 'ALL'])]
    params:
        version = config['Rversion'],
        samplesheet = config['samplesheet']
    log:
        out = 'output/logs/quantNorm.out',
        err = 'output/logs/quantNorm.err'
    shell:
        """
        module load r/{params.version}
        Rscript scripts/quantNorm.R {params.samplesheet} {input} 1> {log.out} 2> {log.err}
        """

rule indexQuant:
    input:
        lambda wildcards: ['output/normquant/{condition}_CPMadjTMM_invNorm.bed'.format(condition=wildcards.condition)]
    output:
        bed = 'output/normquant/{condition}_CPMadjTMM_invNorm.bed.gz',
        index = 'output/normquant/{condition}_CPMadjTMM_invNorm.bed.gz.tbi'
    params:
        version = config['samtoolsVersion']
    log:
        out = 'output/logs/{condition}_indexQuant.out',
        err = 'output/logs/{condition}_indexQuant.err'
    shell:
        """
        module load samtools/{params.version}
        bgzip {input} && tabix -p bed {output.bed} 1> {log.out} 2> {log.err}
        """

rule getPCs:
    input:
        lambda wildcards: ['output/normquant/{condition}_CPMadjTMM_invNorm.bed.gz'.format(condition=wildcards.condition)]
    output:
        'output/covar/{condition}_PC.csv'
    params:
        version = config['Rversion'],
        samplesheet = config['samplesheet']
    log:
        out = 'output/logs/{condition}_getPCs.out',
        err = 'output/logs/{condition}_getPCs.err'
    shell:
        """
        module load r/{params.version}
        Rscript scripts/getPCs.R {params.samplesheet} {input} {wildcards.condition} 1> {log.out} 2> {log.err}
        """

rule makePCcovar:
    input:
        lambda wildcards: ['output/covar/{condition}_PC.csv'.format(condition=wildcards.condition)]
    output:
        'output/covar/{condition}_PCcovar.txt'
    params:
        version = config['Rversion'],
        donorSamplesheet = config['donorSamplesheet']
    log:
        out = 'output/logs/{condition}_makePCcovar.out',
        err = 'output/logs/{condition}_makePCcovar.err'
    shell:
        """
        module load r/{params.version}
        Rscript scripts/formatPCcovariates.R {input} {params.donorSamplesheet} {wildcards.condition} 1> {log.out} 2> {log.err}
        """

# rule PC_eQTL:
#     input:
#         vcf = rules.renameVCFdonors.output.vcf,
#         vcfIndex = rules.renameVCFdonors.output.index,
#         bed = rules.indexQuant.output.bed,
#         bedIndex = rules.indexQuant.output.index,
#         cov = rules.makePCcovar.output
#     output:
#         'output/qtl/PC_{condition}perm1Mb.txt'
#     params:
#         version = config['QTLToolsVersion']
#     log:
#         out = 'output/logs/{condition}_PCeQTL.out',
#         err = 'output/logs/{condition}_PCeQTL.err'
#     shell:
#         """
#         module load qtltools/{params.version}
#         QTLtools cis --vcf {input.vcf} --bed {input.bed} --cov {input.cov} --permute 1000 --window 1000000 --out {output} 1> {log.out} 2> {log.err}
#         """
    
    
        
