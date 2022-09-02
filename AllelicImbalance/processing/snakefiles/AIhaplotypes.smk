#!/usr/bin/env python3
import os
import glob

## Load config file
configfile: "config/config_phasing.yaml"

## Get vcf file path and prefix of VCFproc processed vcf from config file
vcf = config['vcf']
vcf_file = os.path.basename(vcf)
vcf_prefix = vcf_file[:re.search("_nodups_biallelic.vcf.gz", vcf_file).span()[0]]

onsuccess:
    for file in glob.glob('output/AI/chr*_ASEPfinal*_chunk*'):
        os.remove(file)
    for chr in range(1, 23):
        os.remove('genes_chr' + str(chr) + '.done')
        os.remove('output/AI/chr' + str(chr) + '_ASEPheader.txt')
        os.remove('output/AI/chr' + str(chr) + '_geneheader.txt')

rule all:
    input:
        [expand('output/AI/chr{chrom}_ASEPfinal.txt', chrom = range(1, 23))],
        [expand('output/AI/chr{chrom}_ASEPfinal_geneInfo.txt', chrom = range(1, 23))]


rule mergeAlleleCountHaplotypes:
    input:
        alleleCounts = lambda wildcards: 'output/AI/chr{chrom}_alleleCounts_5hets.csv'.format(chrom=wildcards.chrom),
        haps = lambda wildcards: 'output/vcf/' + vcf_prefix + '_nodups_biallelic_chr{chrom}.phased.allele.named.haps'.format(chrom=wildcards.chrom)
    output:
        "output/AI/chr{chrom}_ASEP.csv"
    params:
        chrom = lambda wildcards: wildcards.chrom
    log:
        out = 'output/AI/logs/mergeAlleleCountHaplotypes_chr{chrom}.out'
    shell:
        """
        module load python/3.9.6
        python3 scripts/AIhaplotypes/mergeAlleleCountHaplotypes.py {input.alleleCounts} {input.haps} {params.chrom} 1> {log.out}
        """

rule getHapGenes:
    input:
        rules.mergeAlleleCountHaplotypes.output
    output:
        touch("genes_chr{chrom}.done")
    params:
        chrom = lambda wildcards: wildcards.chrom
    log:
        out = 'output/AI/logs/getHapGenes_chr{chrom}.out'
    shell:
        """
        module load r/4.1.3
        # Get line count for file
        lineNo=`wc -l < {input}`
        Rscript scripts/AIhaplotypes/ASEP_genes.R {input} {params.chrom} ${{lineNo}} 1> {log.out}
        """

rule joinHapASEP:
    input:
        lambda wildcards: 'genes_chr{chrom}.done'.format(chrom=wildcards.chrom)
    output:
        asep = 'output/AI/chr{chrom}_ASEPfinal.txt',
        genes = 'output/AI/chr{chrom}_ASEPfinal_geneInfo.txt'
    params:
        chrom = lambda wildcards: wildcards.chrom
    shell:
        """
        asepFiles=`ls output/AI/chr{params.chrom}_ASEPfinal_chunk*`
        geneFiles=`ls output/AI/chr{params.chrom}_ASEPfinal_geneInfo_chunk*`

        cat output/AI/chr{params.chrom}_ASEPheader.txt ${{asepFiles}} >> {output.asep}
        cat output/AI/chr{params.chrom}_geneheader.txt ${{geneFiles}} >> {output.genes}
        """
