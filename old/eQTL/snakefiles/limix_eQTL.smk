#!/usr/bin/env python3
import pandas as pd
import os, shutil
import re
import glob
import numpy as np
from kneed import KneeLocator

## Load config file
configfile: "config/config_limixQTL_eQTL.yaml"

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


rule_all_inputs = [[expand('output/limix_input/{condition}_featureAnnotation.tsv', condition = ['CTL', 'FNF'])],
                    [expand('output/limix_input/{condition}_phenotypeFile.tsv', condition = ['CTL', 'FNF'])],
                    [expand('output/limix_input/geno.{filext}', filext = ['bed', 'bim', 'fam'])],
                    [expand('output/limix_input/{condition}_PEER_k{Nk}_genoPC_batch_covariateMatrix.tsv', condition = ['CTL', 'FNF'], Nk = n) for n in range(1, 16, 1)],
                    'output/limix_input/sampleMappingFile.tsv',
                    [expand('output/{condition}_limixqtl/qtl_results_all.h5', condition = ['CTL', 'FNF'])],
                    [expand('output/{condition}_limixqtl/snp_metadata_all.h5', condition = ['CTL', 'FNF'])],
                    [expand('output/{condition}_limixqtl/feature_metadata_all.h5', condition = ['CTL', 'FNF'])]]

#----------------------------------------------------------------------
# Rules
#----------------------------------------------------------------------

rule all:
    input:
        rule_all_inputs

include: "eQTL.smk"

rule genoPCA:
    input:
        vcf = 'output/vcf/' + vcf_prefix + '_newcontig_rename.vcf.gz',
        index = 'output/vcf/' + vcf_prefix + '_newcontig_rename.vcf.gz.tbi'
    output:
        pcs = 'output/covar/genotypes.pca',
        pcvar = 'output/covar/genotypes.pca_stats'
    params:
        version = config['QTLToolsVersion']
    log:
        out = 'output/logs/genoPCA.out',
        err = 'output/logs/genoPCA.err'
    shell:
        """
        module load qtltools/{params.version}
        QTLtools pca --vcf {input.vcf} --scale --center --out output/covar/genotypes 1> {log.out} 2> {log.err}
        """

rule genoPCkneedle:
    input:
        rules.genoPCA.output.pcvar
    output:
        'output/covar/genoPCkneedle.txt'
    log:
        out = 'output/logs/genoPCkneedle.out',
        err = 'output/logs/genoPCkneedle.err'
    run:
        PC = pd.read_csv(str(input), header = None, delimiter = " ")

        # Grab var explained
        varExplained = PC.iloc[1].to_numpy()

        # Get rid of 'prop_var' in first index
        varExplained = np.delete(varExplained, 0)

        # Convert to float
        varExplained = varExplained.astype(float)

        # Remove nan
        varExplained_final = varExplained[~np.isnan(varExplained)]

        # Generate array for number of PCs
        x = np.array(range(1, len(PC.columns)))
        x = x[~np.isnan(varExplained)]

        # Calculate knee
        kneedle = KneeLocator(x, varExplained_final, curve = "convex", direction = "decreasing")

        # Write knee to file
        numPCs = kneedle.knee
        file = open(str(output), 'w')
        file.write(str(numPCs))
        file.close()

rule makeFeatureAnnotation:
    input:
        rules.indexQuant.output.bed
    output:
        'output/limix_input/{condition}_featureAnnotation.tsv'
    params:
        version = config['Rversion']
    log:
        out = 'output/logs/{condition}_makeFeatureAnnotation.out',
        err = 'output/logs/{condition}_makeFeatureAnnotation.err'
    shell:
        """
        module load r/{params.version}
        Rscript scripts/limixQTL/makeFeatureAnnotation.R {input} {output} 1> {log.out} 2> {log.err}
        """

rule makePhenotypeFile:
    input:
        rules.indexQuant.output.bed
    output:
        'output/limix_input/{condition}_phenotypeFile.tsv'
    params:
        version = config['Rversion']
    log:
        out = 'output/logs/{condition}_makePhenotypeFile.out',
        err = 'output/logs/{condition}_makePhenotypeFile.err'
    shell:
        """
        module load r/{params.version}
        Rscript scripts/limixQTL/makePhenotypeFile.R {input} {output} 1> {log.out} 2> {log.err}
        """

rule makeGenotypeFile:
    input:
        vcf = rules.renameVCFdonors.output.vcf
    output:
        bed = 'output/limix_input/geno.bed',
        bim = 'output/limix_input/geno.bim',
        fam = 'output/limix_input/geno.fam'
    log:
        out = 'output/logs/makeGenotypeFile.out',
        err = 'output/logs/makeGenotypeFile.err'
    shell:
        """
        module load plink
        plink --vcf {input.vcf} --make-bed --out output/limix_input/geno 1> {log.out} 2> {log.err}
        """

rule makePEERgeno_covariateMatrix:
    input:
        peer = lambda wildcards: 'output/covar/{condition}_PEERfactors_k{Nk}.txt'.format(condition=wildcards.condition, Nk=wildcards.Nk),
        genoPC = rules.genoPCA.output.pcs,
        numgenoPCs = rules.genoPCkneedle.output
    output:
        'output/limix_input/{condition}_PEER_k{Nk}_genoPC_batch_covariateMatrix.tsv'
    log:
        out = 'output/logs/{condition}_makePEERgeno_covariateMatrix_k{Nk}.out',
        err = 'output/logs/{condition}_makePEERgeno_covariateMatrix_k{Nk}.err'
    params:
        version = config['Rversion'],
        donorSamplesheet = config['donorSamplesheet'],
        samplesheet = config['samplesheet']
    shell:
        """
        module load r/{params.version}
        Rscript scripts/limixQTL/makePEERgeno_covariateMatrix.R {input.peer} {params.donorSamplesheet} {params.samplesheet} {wildcards.condition} {wildcards.Nk} {input.genoPC} {input.numgenoPCs} 1> {log.out} 2> {log.err}
        """

rule makeSampleMappingFile:
    input:
        donorSamplesheet = config['donorSamplesheet']
    output:
        'output/limix_input/sampleMappingFile.tsv'
    log:
        out = 'output/logs/makeSampleMappingFile.out',
        err = 'output/logs/makeSampleMappingFile.err'
    params:
        version = config['Rversion']
    shell:
        """
        module load r/{params.version}
        Rscript scripts/limixQTL/makeSampleMappingFile.R {input.donorSamplesheet} {output} 1> {log.out} 2> {log.err}
        """

rule runMainPass:
    input:
        bed = rules.makeGenotypeFile.output.bed,
        bim = rules.makeGenotypeFile.output.bim,
        fam = rules.makeGenotypeFile.output.fam,
        annotation = rules.makeFeatureAnnotation.output
        phenotype = rules.makePhenotypeFile.output
        samples = rules.makeSampleMappingFile.output
    output:
        results = 'output/{condition}_limixqtl/qtl_results_all.h5',
        snp_metadata = 'output/{condition}_limixqtl/snp_metadata_all.h5',
        features = 'output/{condition}_limixqtl/feature_metadata_all.h5'
    log:
        out = 'output/logs/runMainPass.out',
        err = 'output/logs/runMainPass.err'
    params:
        limixPath = config['limixPath'],
        version = config['pythonVersion']
    shell:
        'module load python/{params.version} &&'
        'python3 {params.limixPath}/run_QTL_analysis.py '
        '--plink output/limix_input/geno '
        '--annotation_file {input.annotation} '
        '--phenotype_file {input.phenotype} '
        '--window 1000000 '
        '--genomic_range all '
        '--number_of_permutations 1000 '
        '--sample_mapping_file {input.samples} '
        '--cis '
        '--output_directory output/{condition}_limixqtl 1> {log.out} 2> {log.err}'

#rule InteractionPass: