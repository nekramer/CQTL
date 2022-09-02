#!/usr/bin/env python3

rule all:
    input:
        'data/ASEP_results.csv'

rule runASEP:
    input:
        lambda wildcards: '../processing/output/AI/chr{chrom}_ASEPfinal.txt'.format(chrom=wildcards.chrom)
    output:
        temp('data/chr{chrom}_ASEP.rda')
    params:
        chrom = lambda wildcards: wildcards.chrom
    threads: 8
    log:
        out = "logs/runASEP_chr{chrom}.out"
    shell:
        """
        module load r/4.1.3
        Rscript scripts/runASEP.R {input} {params.chrom} 1> {log.out}
        """
    
rule combineASEP:
    input:
        [expand('data/chr{chrom}_ASEP.rda', chrom = range(1, 23))]
    output:
        'data/ASEP_results.csv'
    log:
        out = "logs/combineASEP.out"
    shell:
        """
        module load r/4.1.3
        Rscript scripts/combineASEP.R
        """