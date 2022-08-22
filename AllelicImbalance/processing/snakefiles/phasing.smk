#!/usr/bin/env python3

## Load config file
configfile: "config/config_phasing.yaml"

## Get vcf file path and prefix of VCFproc processed vcf from config file
vcf = config['vcf']
vcf_file = os.path.basename(vcf)
vcf_prefix = vcf_file[:re.search("_nodups_biallelic.vcf.gz", vcf_file).span()[0]]

rule all:
    expand('output/vcf/' + vcf_prefix + '_nodups_biallelic_chr{chr}.phased.haps', chr = range(1, 23)),
    expand('output/vcf/' + vcf_prefix + '_nodups_biallelic_chr{chr}.phased.sample', chr = range(1, 23))

rule splitVCF:
    input: 
        'output/vcf/' + vcf_prefix + '_nodups_biallelic.vcf.gz'
    output:
        expand('output/vcf/' + vcf_prefix  + '_nodups_biallelic_chr{chr}', chr = range(1, 23))
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
        lambda wildcards: 'output/vcf/' + vcf_prefix + '_nodups_biallelic_chr{chr}.vcf.gz'.format(chr=wildcards.chr)
    output:
       haps = 'output/vcf/' + vcf_prefix + '_nodups_biallelic_chr{chr}.phased.haps',
       samples = 'output/vcf/' + vcf_prefix + '_nodups_biallelic_chr{chr}.phased.sample'
    params:
        geneticMapDir = config['geneticMapDir'],
        geneticMapPrefix = config['geneticMapPrefix'],
        chrom = lambda wildcards: wildcards.chr
    threads: 8
    log:
        out = 'output/vcf/logs/phaseData_chr{chr}.out'
    shell:
        """
        module load shapeit
        shapeit --input-vcf {input} --input-map {params.geneticMapDir}/{params.geneticMapPrefix}{params.chrom}_sorted.txt --thread {threads} --output {output.haps} {output.samples} --output-log {log.out}
        """

