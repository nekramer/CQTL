#!/usr/bin/env python3

## Load config file
configfile: "config/config_phasing.yaml"

## Get vcf file path and prefix of VCFproc processed vcf from config file
vcf = config['vcf']
vcf_file = os.path.basename(vcf)
vcf_prefix = vcf_file[:re.search("_nodups_biallelic.vcf.gz", vcf_file).span()[0]]

rule all:
    input:
        [expand('output/vcf/{prefix}_nodups_biallelic_chr{chrom}.phased.{ext}', prefix = vcf_prefix, chrom = range(1, 23), ext = ['haps', 'sample'])]

rule splitVCF:
    input: 
        'output/vcf/' + vcf_prefix + '_nodups_biallelic.vcf.gz'
    output:
        expand('output/vcf/' + vcf_prefix  + '_nodups_biallelic_chr{chrom}', chrom = range(1, 23))
    params: 
        prefix = "output/vcf/" + vcf_prefix + "_nodups_biallelic_chr"
    log:
        out = 'output/vcf/logs/splitVCF.out'
    shell:
        """
        module load samtools
        for chr in {1..22}
        do
            bcftools view {input} --regions chr${chr} -o {params.prefix}${chr} -Oz 1> {log.out}
        done
        """

rule phaseData:
    input:
        lambda wildcards: 'output/vcf/' + vcf_prefix + '_nodups_biallelic_chr{chrom}.vcf.gz'.format(chrom=wildcards.chrom)
    output:
       haps = 'output/vcf/' + vcf_prefix + '_nodups_biallelic_chr{chrom}.phased.haps',
       samples = 'output/vcf/' + vcf_prefix + '_nodups_biallelic_chr{chrom}.phased.sample'
    params:
        geneticMapDir = config['geneticMapDir'],
        geneticMapPrefix = config['geneticMapPrefix'],
        chrom = lambda wildcards: wildcards.chrom
    threads: 8
    log:
        out = 'output/vcf/logs/phaseData_chr{chrom}.out'
    shell:
        """
        module load shapeit
        shapeit --input-vcf {input} --input-map {params.geneticMapDir}/{params.geneticMapPrefix}{params.chrom}_sorted.txt --thread {threads} --output-max {output.haps} {output.samples} --output-log {log.out}
        """

