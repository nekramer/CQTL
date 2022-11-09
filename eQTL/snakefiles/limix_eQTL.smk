#!/usr/bin/env python3
import pandas as pd
import os, shutil
import re
import glob
import numpy as np
from kneed import KneeLocator

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


rule_all_inputs = []



#----------------------------------------------------------------------
# Rules
#----------------------------------------------------------------------

rule all:
    input:
        rule_all_inputs

include: "../../rules/VCFprocessing.smk"

include: "../../rules/RNAprocessing.smk"

rule renameVCFdonors:
    input:
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

rule getPEER:
    input:
        lambda wildcards: ['output/normquant/{condition}_CPMadjTMM_invNorm.bed.gz'.format(condition=wildcards.condition)]
    output:
        [expand('output/covar/{{condition}}_PEERfactors_k{Nk}.txt', Nk = n) for n in range(1, 16, 1)]
    log:
        out = 'output/logs/{condition}_getPEER.out',
        err = 'output/logs/{condition}_getPEER.err'
    shell:
        """
        module load r/4.2.0
        Rscript scripts/PEERfactors.R {input} {wildcards.condition} 1> {log.out} 2> {log.err}
        """

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

