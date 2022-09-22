#!/usr/bin/env python3

import pandas as pd
import os, re
import shutil

## Load config file
configfile: "config/config_AIprocessing.yaml"

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

## Extract grouped donors, getting unique values for each key
groupedDonors = samples.groupby('Donor')['Sample'].apply(list).to_dict()
for key in groupedDonors.keys():
    groupedDonors[key] = tuple(set(groupedDonors[key]))

## Get vcf file path and prefix of VCFproc processed vcf from config file
vcf = config['vcf']
vcf_file = os.path.basename(vcf)
vcf_prefix = vcf_file[:re.search("_ALL_qc.vcf.gz", vcf_file).span()[0]]

onsuccess:
    # Remove intermediate sample joined allele count files
    for key in read1:
        os.remove('output/{group}/alleleCounts/{group}_alleleCounts_joined.csv'.format(group = key))

    # Remove intermediate joined and checked allele counts files and genotyping counts for each donor
    for key in groupedDonors:
        shutil.rmtree('output/AI/{donor}'.format(donor = key))

    # Remove additional subsetted vcf files
    os.remove('output/vcf/' + vcf_prefix + '_nodups_biallelic_AI.recode.vcf.gz')
    os.remove('output/vcf/' + vcf_prefix + '_nodups_biallelic_AI.recode.vcf.gz.tbi')

rule all:
    input:
        [expand("output/{group}/alleleCounts/{group}_alleleCounts_joined.csv", group = key) for key in read1],
        "output/vcf/" + vcf_prefix + "_nodups_biallelic_AI.recode.vcf.gz",
        "output/vcf/" + vcf_prefix + "_nodups_biallelic_AI.recode.vcf.gz.tbi",
        [expand("output/AI/{donor}/{donor}_alleleCounts_joined.csv", donor = key) for key in groupedDonors],
        [expand('output/AI/{donor}/{donor}_genoCounts.csv', donor = key) for key in groupedDonors],
        [expand('output/AI/{donor}/{donor}_alleleCounts_checked.csv', donor = key) for key in groupedDonors],
        'output/AI/alleleCounts.csv',
        'output/AI/numVariantHets.csv',
        'output/AI/alleleCountsplits.txt'

include: "../../../rules/VCFprocessing.smk"

include: "../../../rules/RNAprocessing.smk"

rule assignGroups:
    input:
        R = rules.align.output,
        I = rules.index.output
    output:
        "output/{group}/grouped/{group}.grouped.sort.bam"
    threads: 1
    log:
        err = "output/{group}/logs/{group}_assignGroups.err"
    shell:
        """
        mkdir -p output/{wildcards.group}/grouped
        module load picard/2.2.4
        module load java/10.0.2
        java -jar /nas/longleaf/apps/picard/2.2.4/picard-tools-2.2.4/picard.jar AddOrReplaceReadGroups I={input.R} O={output} RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM={wildcards.group} SORT_ORDER=coordinate 2> {log.err}
        """

rule countReads:
    input:
        bam = rules.assignGroups.output,
        vcf = rules.zipVCF2.output,
        i = rules.indexVCF2.output
    output:
        "output/{group}/alleleCounts/{group}_alleleCounts.csv"
    threads: 1
    log:
        err = "output/{group}/logs/{group}_countReads.err"
    params:
        sequence = config['sequence']
    shell:
        """
        mkdir -p output/{wildcards.group}/alleleCounts
        module load gatk/4.1.7.0
        module load python/3.6.6
        gatk ASEReadCounter --input {input.bam} --variant {input.vcf} --output {output} --output-format RTABLE --min-base-quality 20 --disable-read-filter NotDuplicateReadFilter --reference {params.sequence} 2> {log.err}
        """

rule editDonors:
    input:
        vcf = rules.zipVCF2.output
    output:
        temp(touch('editDonors.done'))
    params:
        donors = ",".join(samples['Donor'].unique().tolist())
    shell:
        """
        module load samtools
        module load python/3.9.6
        bcftools query -l {input.vcf} > donors.txt

        python3 scripts/RNAproc/matchDonors.py donors.txt {params.donors}
        """

rule getVariants:
    input: 
        [expand("output/{group}/alleleCounts/{group}_alleleCounts.csv", group = key) for key in read1]
    output:
        temp('output/AI/variants.csv')
    log:
        out = 'output/AI/logs/getVariants.out'
    shell:
        """
        module load r/4.2.1
        mkdir -p output/AI/logs
        Rscript scripts/AIprocessing/getVariants.r {input} 1> {log.out}
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
        python3 scripts/AIprocessing/mergeSampleVariants.py {input} {params.minTotalAlleleCounts} {params.minAlleleCounts} 1> {log.out}
        """

rule concatDonorConditions:
    input:
        lambda wildcards: expand('output/{group}/alleleCounts/{group}_alleleCounts_joined.csv', group = groupedDonors[wildcards.donor])
    output:
        'output/AI/{donor}/{donor}_alleleCounts_joined.csv'
    params:
        donor = lambda wildcards: wildcards.donor
    log:
        out = 'output/AI/logs/{donor}_concatDonorConditions.out'
    shell:
        """
        module load python/3.9.6
        mkdir -p output/AI/logs
        mkdir -p output/AI/{params.donor}
        python3 scripts/AIprocessing/concatDonorConditions.py {params.donor} {input} 1> {log.out}
        """
        
rule VCFoverlapVariants:
    input:
        vcf = rules.zipVCF2.output,
        variants = "output/AI/variants.csv"
    output:
        "output/vcf/" + vcf_prefix + "_nodups_biallelic_AI.recode.vcf.gz"
    params:
        prefix = "output/vcf/" + vcf_prefix + "_nodups_biallelic_AI"
    shell:
        """
        module load vcftools
        module load samtools
        vcftools --gzvcf {input.vcf} --snps {input.variants} --recode --recode-INFO-all --out {params.prefix}
        bgzip {params.prefix}.recode.vcf   
        """

rule VCFoverlapVariantsIndex:
    input:
        rules.VCFoverlapVariants.output
    output:
        "output/vcf/" + vcf_prefix + "_nodups_biallelic_AI.recode.vcf.gz.tbi"
    log:
        out = "output/vcf/logs/VCFoverlapVariantsIndex.out",
        err = "output/vcf/logs/VCFoverlapVariantsIndex.err"
    shell:
        """
        module load samtools
        tabix -p vcf {input} 2> {log.err} 1> {log.out}
        """

rule getGenoCounts:
    input:
        vcf = rules.VCFoverlapVariants.output,
        index = rules.VCFoverlapVariantsIndex.output,
        donorConversions = 'donors.txt',
        donorDone = 'editDonors.done'
    output:
        'output/AI/{donor}/{donor}_genoCounts.csv'
    params:
        donor = lambda wildcards: wildcards.donor
    log:
        out = 'output/AI/logs/{donor}_getGenoCounts.out'
    shell:
        """
        module load python/3.9.6
        mkdir -p output/AI/logs
        mkdir -p output/AI/{params.donor}
        python3 scripts/AIprocessing/getGenoCounts.py {input.vcf} {params.donor} {input.donorConversions} 1> {log.out}
        """

rule checkDonorVariants:
    input: 
        lambda wildcards: ['output/AI/{donor}/{donor}_alleleCounts_joined.csv'.format(donor=wildcards.donor)],
        lambda wildcards: ['output/AI/{donor}/{donor}_genoCounts.csv'.format(donor=wildcards.donor)]
    output:
        'output/AI/{donor}/{donor}_alleleCounts_checked.csv'
    params:
        donor = lambda wildcards: wildcards.donor
    log:
        out = 'output/AI/logs/{donor}_checkDonorVariants.out'
    shell:
        """
        module load python/3.9.6
        mkdir -p output/AI/logs
        python3 scripts/AIprocessing/checkDonorVariants.py {input} {params.donor} 1> {log.out}
        """

rule concatAlleleCounts:
    input: 
        [expand("output/AI/{donor}/{donor}_alleleCounts_checked.csv", donor = key) for key in groupedDonors]
    output:
        'output/AI/alleleCounts.csv'
    log:
        out = 'output/AI/logs/concatAlleleCounts.out'
    shell:
        """
        module load python/3.9.6
        mkdir -p output/AI/logs
        python3 scripts/AIprocessing/concatAlleleCounts.py {input} 1> {log.out}
        """

rule checkVariantHets:
    input:
        rules.concatAlleleCounts.output
    output:
        'output/AI/numVariantHets.csv'
    log:
        out = 'output/AI/logs/checkVariantHets.out'
    shell:
        """
        module load python/3.9.6
        mkdir -p output/AI/logs
        python3 scripts/AIprocessing/checkVariantHets.py {input} 1> {log.out}
        """

rule splitAlleleCountsFile:
    input:
        rules.concatAlleleCounts.output
    output:
        'output/AI/alleleCountsplits.txt'
    log:
        err = 'output/AI/logs/splitAlleleCountsFile.err'
    shell:
        """
        mkdir -p output/AI/alleleCountSplits
        split {input} -l 100000 output/AI/alleleCountSplits/alleleCounts_split 2> {log.err}
        ls -1 output/AI/alleleCountSplits > {output} 2>> {log.err}
        """

