#!/usr/bin/env python3

## Load config file
configfile: "config/config_phasing.yaml"

## Get vcf file path and prefix of VCFproc processed vcf from config file
vcf = config['vcf']
vcf_file = os.path.basename(vcf)
vcf_prefix = vcf_file[:re.search("_nodups_biallelic.vcf.gz", vcf_file).span()[0]]

rule all:
    input:
        [expand('output/vcf/{prefix}_nodups_biallelic_chr{chrom}.phased.{ext}', prefix = vcf_prefix, chrom = range(1, 23), ext = ['haps', 'sample'])],
        [expand('output/vcf/{prefix}_nodups_biallelic_chr{chrom}.phased.allele.named.haps', prefix = vcf_prefix, chrom = range(1, 23))]

rule splitVCF:
    input: 
        'output/vcf/' + vcf_prefix + '_nodups_biallelic.vcf.gz'
    output:
        expand('output/vcf/' + vcf_prefix  + '_nodups_biallelic_chr{chrom}.vcf.gz', chrom = range(1, 23))
    params: 
        prefix = "output/vcf/" + vcf_prefix + "_nodups_biallelic_chr"
    log:
        out = 'output/vcf/logs/splitVCF.out'
    shell:
        """
        module load samtools
        for chr in {{1..22}}
        do
            bcftools view {input} --regions chr${{chr}} -o {params.prefix}${{chr}} -Oz 1> {log.out}
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
        module load shapeit/2.837
        shapeit --input-vcf {input} --input-map {params.geneticMapDir}/{params.geneticMapPrefix}{params.chrom}_sorted.txt --thread {threads} --output-max {output.haps} {output.samples} --output-log {log.out}
        """

rule hapsToAlleles:
    input: 
        rules.phaseData.output.haps
    output:
        temp('output/vcf/' + vcf_prefix + '_nodups_biallelic_chr{chrom}.phased.allele.haps')
    params:
        chrom = lambda wildcards: wildcards.chrom,
        prefix = vcf_prefix
    log:
        out = 'output/vcf/logs/haptToAlleles_chr{chrom}.out'
    shell:
        """
        module load python/3.9.6
        python3 scripts/phasing/hapsToAlleles.py {input} {params.chrom} {params.prefix} 1> {log.out}
        """

rule renameHapColumns:
    input:
        haps = rules.hapsToAlleles.output,
        samples = rules.phaseData.output.samples
    output:
        'output/vcf/' + vcf_prefix + '_nodups_biallelic_chr{chrom}.phased.allele.named.haps'
    params:
        chrom = lambda wildcards: wildcards.chrom,
        prefix = vcf_prefix
    log:
        out = 'output/vcf/logs/renameHapColumns_chr{chrom}.out'
    shell:
        """
        module load python/3.9.6
        python3 scripts/phasing/renameHapColumns.py {input.haps} {input.samples} {params.chrom} {params.prefix} 1> {log.out}
        """
