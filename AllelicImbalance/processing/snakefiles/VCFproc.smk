#!/usr/bin/env python3
import os, re

## Load config file
configfile: "config/config.yaml"

## Get file path of post-imputed, qc'd gzipped vcf file
vcf = config["vcf"]

## Determine a prefix for the vcf file
vcf_file = os.path.basename(vcf)
vcf_prefix = vcf_file[:re.search("_ALL_qc.vcf.gz", vcf_file).span()[0]]

## Define actions on success
onsuccess:

    ## Success message
    print("VCFs processed successfully for Allelic Imbalance!")

    ## Remove dummy output file
    os.remove("edit.done")


rule all:
    input:
        'output/vcf/' + vcf_prefix + '_nodups_biallelic.vcf.gz',
        'output/vcf/' + vcf_prefix + '_nodups_biallelic.vcf.gz.tbi',
        'edit.done'

include: "../../../rules/VCFprocessing.smk"

rule selectVariants:
    input:
        v = rules.updateConfig.output.v,
        i = rules.updateConfig.output.i
    output:
        v = temp('output/vcf/' + vcf_prefix + '_biallelic.vcf.gz'),
        i = temp('output/vcf/' + vcf_prefix + '_biallelic.vcf.gz.tbi')
    threads: 1
    log:
        out = "output/vcf/logs/selectVariants.out",
        err = "output/vcf/logs/selectVariants.err"
    params:
        sequence = config['sequence']
    shell:
        """
        module load gatk/4.1.7.0
        gatk SelectVariants --variant {input.v} -R {params.sequence} --select-type-to-include SNP -O {output.v} --restrict-alleles-to BIALLELIC 2> {log.err} 1> {log.out}
        """

rule removeDuplicates:
    input:
        v = rules.selectVariants.output.v,
        i = rules.selectVariants.output.i
    output:
        'output/vcf/' + vcf_prefix + '_nodups_biallelic.vcf'
    threads: 4
    log:
        err = "output/vcf/logs/removeDuplicates.err"
    shell:
        """
        module load samtools
        bcftools norm -d any --threads {threads} -o {output} {input.v} 2> {log.err}
        """

rule zipVCF2:
    input:
        rules.removeDuplicates.output
    output:
        'output/vcf/' + vcf_prefix + '_nodups_biallelic.vcf.gz'
    threads: 4
    log:
        err = "output/vcf/logs/zipVCF2.err"
    shell:
        """
        module load samtools
        bgzip --threads {threads} {input} 2> {log.err}
        """  

rule indexVCF2:
    input:
        rules.zipVCF2.output
    output:
        'output/vcf/' + vcf_prefix + '_nodups_biallelic.vcf.gz.tbi'
    threads: 1
    log:
        out = "output/vcf/logs/indexVCF2.out",
        err = "output/vcf/logs/indexVCF2.err"
    shell:
        """
        module load samtools
        tabix -p vcf {input} 2> {log.err} 1> {log.out}
        """

rule addVCFPath:
    input:
        name = rules.zipVCF2.output,
        config = 'config/config_AIanalysis.yaml'
    output:
        touch('edit.done')
    shell:
        """
        sed -i 's#vcf:#& "'{input.name}'"#' {input.config}
        """