#!/usr/bin/env python3

## Load config file
configfile: "config/config_reQTL.yaml"
vcf = config['vcf']
rna = config['rna']
Nk = config['PEERfactors']

batches = ['RNAKitBatch', 'RNASequencingBatch', 'genoBatch', 'DNAKitBatch']
peerCov = 'output/covar/ALL_PEERk{Nk}_genoPC'
fileExt = ''

for b in batches:
    b_include = config[b]
    if b_include == "TRUE":
        peerCov += '_{}'.format(b)
        fileExt += '_{}'.format(b)

peerCov += ".csv"

rule all:
    input:
        expand('output/reQTL/{condition}_sig_reQTLs.rds', condition = ['CTL', 'FNF']),
        expand('output/reQTL/{condition}_sig_reQTLs.csv', condition = ['CTL', 'FNF'])

# Separate eGenes that are only found in CTL vs only found in FNF
rule sep_eGenes:
    input:
        CTL_eQTL = config['eQTL'] + 'CTL_PEER_k' + str(Nk) + '_genoPC' + fileExt + '_perm1Mb_FDR.txt',
        FNF_eQTL = config['eQTL'] + 'FNF_PEER_k' + str(Nk) + '_genoPC' + fileExt + '_perm1Mb_FDR.txt'   
    output:
        'output/reQTL/CTLonly_sig_eGenes.csv',
        'output/reQTL/FNFonly_sig_eGenes.csv'
    params:
        version = config['Rversion'],
        correction = config['correction'],
        threshold = config['FDRthreshold']
    log:
        out = 'output/logs/sep_eGenes.out',
        err = 'output/logs/sep_eGenes.err'
    shell:
        """
        module load r/{params.version}
        Rscript scripts/responseQTL/separate_eGenes.R {input.CTL_eQTL} {input.FNF_eQTL} {params.correction} {params.threshold} 1> {log.out} 2> {log.err}
        """ 

# Filter the VCF file for both sets of significant lead variants
rule subsetVCF_leadvar:
    input:
        eGenes_CTL = 'output/reQTL/CTLonly_sig_eGenes.csv',
        eGenes_FNF = 'output/reQTL/FNFonly_sig_eGenes.csv',
        vcf = vcf
    output:
        'output/reQTL/ALL_leadVars.vcf.gz'
    params:
        version = config['gatkVersion']
    log:
        out = 'output/logs/subsetVCF_leadvar.out',
        err = 'output/logs/subsetVCF_leadvar.err'
    shell:
        """
        module load gatk/{params.version}
        # Grab variantID from each condition-specific eGenes list, removing header
        cut -d "," -f7 {input.eGenes_CTL} > output/reQTL/CTL_variants.list
        sed -i '1d' output/reQTL/CTL_variants.list
        cut -d "," -f7 {input.eGenes_FNF} > output/reQTL/FNF_variants.list
        sed -i '1d' output/reQTL/FNF_variants.list
        # Concatenate
        cat output/reQTL/CTL_variants.list output/reQTL/FNF_variants.list > output/reQTL/ALL_variants.list
        # Subset vcf for these
        gatk SelectVariants -V {input.vcf} --keep-ids output/reQTL/ALL_variants.list -O {output} 1> {log.out} 2> {log.err}
        """

# Get PEER factors for ALL norm RNA counts
rule getPEER_ALL:
    input:
        rna
    output:
        factors = [expand('output/covar/ALL_PEERfactors_k{Nk}.txt', Nk = Nk)],
        var = [expand('output/covar/ALL_PEERfactors_k{Nk}_variance.txt', Nk = Nk)]
    params:
        Nk = config['PEERfactors']
    log:
        out = 'output/logs/ALL_getPEER.out',
        err = 'output/logs/ALL_getPEER.err'
    shell:
        """
        module load r/4.2.0
        Rscript scripts/PEERfactors.R {input} ALL {params.Nk} FALSE 1> {log.out} 2> {log.err}
        """

# Organize covariates for ALL RNA counts
rule makePEERcovar_geno_ALL:
    input:
        peer = rules.getPEER_ALL.output.factors
    output:
        [expand(peerCov, Nk = Nk)]
    params:
        version = config['Rversion'],
        genoPC = config['genoPC'],
        genoPCkneedle = config['genoPCkneedle'],
        donorSamplesheet = config['donorSamplesheet'],
        dnaSamplesheet = config['dnaSamplesheet'],
        samplesheet = config['samplesheet'],
        RNAKitBatch = config['RNAKitBatch'],
        RNASequencingBatch = config['RNASequencingBatch'],
        genoBatch = config['genoBatch'],
        DNAKitBatch = config['DNAKitBatch'],
        Nk = Nk
    log:
        out = 'output/logs/makePEERcovar_geno_ALL.out',
        err = 'output/logs/makePEERcovar_geno_ALL.err'
    shell:
        """
        module load r/{params.version}
        numgenoPCs=`cat {params.genoPCkneedle}`

        Rscript scripts/responseQTL/formatPEERcovariates_geno_ALL.R {input.peer} {params.donorSamplesheet} {params.dnaSamplesheet} {params.genoBatch} {params.DNAKitBatch} {params.samplesheet} {params.RNAKitBatch} {params.RNASequencingBatch} {params.genoPC} ${{numgenoPCs}} {output} 1> {log.out} 2> {log.err}
        """

rule get_reQTLs:
    input:
        normExpression = rna,
        covariates = rules.makePEERcovar_geno_ALL.output,
        vcf = rules.subsetVCF_leadvar.output,
        eGene = lambda wildcards: 'output/reQTL/{condition}only_sig_eGenes.csv'.format(condition = wildcards.condition)
    output:
        rds = 'output/reQTL/{condition}_sig_reQTLs.rds',
        csv = 'output/reQTL/{condition}_sig_reQTLs.csv'
    params:
        version = config['Rversion'],
        threshold = config['FDRthreshold']
    log:
        out = "output/logs/get{condition}reQTLs.out",
        err = "output/logs/get{condition}reQTLs.err"
    shell:
        """
        module load r/{params.version}
        Rscript scripts/responseQTL/testInteraction.R {input.normExpression} {input.covariates} {input.vcf} {input.eGene} {params.threshold} {output.rds} {output.csv} 1> {log.out} 2> {log.err}
        """

# Get rsids

# plot side-by-side
# plot individually