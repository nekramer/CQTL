#!/usr/bin/env python3
import pandas as pd
import os, shutil
import re
import glob
import numpy as np
from kneed import KneeLocator


## Load config file
configfile: "config/config_delta_eQTL.yaml"

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
        "output/normquant/ALL_CPMadjTMM_invNorm.bed"
    params:
        version = config['Rversion'],
        samplesheet = config['samplesheet']
    log:
        out = 'output/logs/quantNorm.out',
        err = 'output/logs/quantNorm.err'
    shell:
        """
        module load r/{params.version}
        Rscript scripts/quantNorm.R {params.samplesheet} ALL {input} 1> {log.out} 2> {log.err}
        """

rule getPEER:
    input:
        'output/normquant/ALL_CPMadjTMM_invNorm.bed.gz'
    output:
        [expand('output/covar/ALL_PEERfactors_k{Nk}.txt', Nk = n) for n in range(1, 16, 1)]
    log:
        out = 'output/logs/delta_eQTL/getPEER.out',
        err = 'output/logs/delta_eQTL/getPEER.err'
    shell:
        """
        module load r/4.2.0
        Rscript scripts/PEERfactors.R {input} ALL 1> {log.out} 2> {log.err}
        """

rule checkPEER:
    input:
        [expand("output/covar/ALL_PEERfactors_k{Nk}.txt", Nk = n) for n in range(1, 16, 1)]
    output:
        'output/covar/ALL_validPEER.txt'
    log:
        out = 'output/logs/ALL_checkPEER.out',
        err = 'output/logs/ALL_checkPEER.err'
    run:
        for file in input:
            peerFactors = pd.read_csv(str(file))
            peerFactors = peerFactors.drop('Donor', axis = 1)
            if not peerFactors.isnull().values.all():
                index = input.index(file)
                Nk = range(1, 16, 1)
                # Write Nk as a valid PEER factor
                outFile = open(str(output), 'a')
                outFile.write(str(Nk[index])+"\n")
                outFile.close()

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


rule formatCovariates:
    input:
        peer = lambda wildcards: 'output/covar/ALL_PEERfactors_k{Nk}.txt'.format(Nk=wildcards.Nk),
        check = 'output/covar/ALL_validPEER.txt',
        genoPC = rules.genoPCA.output.pcs,
        numgenoPCs = rules.genoPCkneedle.output
    output:
        'output/covar/ALL_PEER_k{Nk}_genoPC_batch_covar.txt'
    params:
        version = config['Rversion'],
        donorSamplesheet = config['donorSamplesheet'],
        samplesheet = config['samplesheet']
    log:
        out = 'output/logs/ALL_makePEERcovar_geno_k{Nk}.out',
        err = 'output/logs/ALL_makePEERcovar_geno_k{Nk}.err'
    shell:
        """
        module load r/{params.version}
        numgenoPCs=`cat {input.numgenoPCs}`
        Rscript scripts/formatPEERcovariates_geno.R {input.peer} {params.donorSamplesheet} {params.samplesheet} {wildcards.condition} {wildcards.Nk} {input.check} TRUE {input.genoPC} ${{numgenoPCs}} 1> {log.out} 2> {log.err}
        """

