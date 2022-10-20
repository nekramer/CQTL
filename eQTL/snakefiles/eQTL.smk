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

rule_all_inputs = [[expand("output/covar/{condition}_PC.csv", condition = ['CTL', 'FNF', 'ALL'])],
        [expand('output/pca/{condition}_PC.csv', condition = ['CTL', 'FNF'])],
        [expand('output/covar/{condition}_PCkneedle.txt', condition = ['CTL', 'FNF'])],
        [expand('output/covar/{condition}_PEERfactors_k{Nk}.txt', condition = ['CTL', 'FNF'], Nk = n) for n in range(10, 81, 10)],
        [expand('output/plots/{condition}_PEER{Nk}_correlation.png', condition = ['CTL', 'FNF'], Nk = n) for n in range(10, 81, 10)],
        [expand('output/covar/{condition}_validPEER.txt', condition = ['CTL', 'FNF'])],
        [expand('output/plots/{condition}_pca.pdf', condition = ['CTL', 'FNF', 'ALL'])],
        [expand('output/plots/{condition}_screeplot.pdf', condition = ['CTL', 'FNF', 'ALL'])],
        [expand('output/plots/{condition}_PC_correlation.png', condition = ['CTL', 'FNF', 'ALL'])]]

#[expand("output/mbv/{group}.bamstat.txt", group = key) for key in read1]

if config['genoCovar'] == 'yes':
    include: 'genoCovariate.smk'
    if config['batchCovar'] == 'TRUE':
        rule_all_inputs.extend([[expand('output/covar/{condition}_PC_genoPC_batch_covar.txt', condition = ['CTL', 'FNF'])]])
        rule_all_inputs.extend([[expand('output/covar/{condition}_PEER_k{Nk}_genoPC_batch_covar.txt', condition = ['CTL', 'FNF'], Nk = n) for n in range(10, 81, 10)]])
        pcCov = ['output/covar/{condition}_PC_genoPC_batch_covar.txt']
        peerCov = ['output/covar/{condition}_PEER_k{Nk}_genoPC_batch_covar.txt']
        rule_all_inputs.extend([[expand('output/qtl/{condition}_PC_genoPC_batch_perm1Mb.txt', condition = ['CTL', 'FNF'])]])
        rule_all_inputs.extend([[expand('output/qtl/{condition}_PEER_k{Nk}_genoPC_batch_perm1Mb.txt', condition = ['CTL', 'FNF'], Nk = n) for n in range(10, 81, 10)]])
        pcQTL = ['output/qtl/{condition}_PC_genoPC_batch_perm1Mb.txt']
        peerQTL = ['output/qtl/{condition}_PEER_k{Nk}_genoPC_batch_perm1Mb.txt']
        rule_all_inputs.extend([[expand('output/qtl/{condition}_PC_genoPC_batch_perm1Mb_FDR.txt', condition = ['CTL', 'FNF'])]])
        rule_all_inputs.extend([[expand('output/qtl/{condition}_PEER_k{Nk}_genoPC_batch_perm1Mb_FDR.txt', condition = ['CTL', 'FNF'], Nk = n) for n in range(10, 81, 10)]])
        pcMultipleTesting = ['output/qtl/{condition}_PC_genoPC_batch_perm1Mb_FDR.txt']
        peerMultipleTesting = ['output/qtl/{condition}_PEER_k{Nk}_genoPC_batch_perm1Mb_FDR.txt']

        peer_eGene = [expand('output/qtl/{{condition}}_PEER_k{Nk}_genoPC_batch_perm1Mb_FDR.txt', Nk = n) for n in range(10, 81, 10)]

        rule_all_inputs.extend([[expand('output/summary/{condition}_genoPC_batch.txt', condition = ['CTL', 'FNF'])]])
        eGene = ['output/summary/{condition}_genoPC_batch.txt']
    else:
        rule_all_inputs.extend([[expand('output/covar/{condition}_PC_genoPC_covar.txt', condition = ['CTL', 'FNF'])]])
        rule_all_inputs.extend([[expand('output/covar/{condition}_PEER_k{Nk}_genoPC_covar.txt', condition = ['CTL', 'FNF'], Nk = n) for n in range(10, 81, 10)]])
        pcCov = ['output/covar/{condition}_PC_genoPC_covar.txt']
        peerCov = ['output/covar/{condition}_PEER_k{Nk}_genoPC_covar.txt']
        rule_all_inputs.extend([[expand('output/qtl/{condition}_PC_genoPC_perm1Mb.txt', condition = ['CTL', 'FNF'])]])
        rule_all_inputs.extend([[expand('output/qtl/{condition}_PEER_k{Nk}_genoPC_perm1Mb.txt', condition = ['CTL', 'FNF'], Nk = n) for n in range(10, 81, 10)]])
        pcQTL = ['output/qtl/{condition}_PC_genoPC_perm1Mb.txt']
        peerQTL = ['output/qtl/{condition}_PEER_k{Nk}_genoPC_perm1Mb.txt']
        rule_all_inputs.extend([[expand('output/qtl/{condition}_PC_genoPC_perm1Mb_FDR.txt', condition = ['CTL', 'FNF'])]])
        rule_all_inputs.extend([[expand('output/qtl/{condition}_PEER_k{Nk}_genoPC_perm1Mb_FDR.txt', condition = ['CTL', 'FNF'], Nk = n) for n in range(10, 81, 10)]])
        pcMultipleTesting = ['output/qtl/{condition}_PC_genoPC_perm1Mb_FDR.txt']
        peerMultipleTesting = ['output/qtl/{condition}_PEER_k{Nk}_genoPC_perm1Mb_FDR.txt']

        peer_eGene = [expand('output/qtl/{{condition}}_PEER_k{Nk}_genoPC_perm1Mb_FDR.txt', Nk = n) for n in range(10, 81, 10)]

        rule_all_inputs.extend([[expand('output/summary/{condition}_genoPC.txt', condition = ['CTL', 'FNF'])]])
        eGene = ['output/summary/{condition}_genoPC.txt']
else:
    include: 'no_genoCovariate.smk'
    if config['batchCovar'] == 'TRUE':
        rule_all_inputs.extend([[expand('output/covar/{condition}_PC_batch_covar.txt', condition = ['CTL', 'FNF'])]])
        rule_all_inputs.extend([[expand('output/covar/{condition}_PEER_k{Nk}_batch_covar.txt', condition = ['CTL', 'FNF'], Nk = n) for n in range(10, 81, 10)]])
        pcCov = ['output/covar/{condition}_PC_batch_covar.txt']
        peerCov = ['output/covar/{condition}_PEER_k{Nk}_batch_covar.txt']
        rule_all_inputs.extend([[expand('output/qtl/{condition}_PC_batch_perm1Mb.txt', condition = ['CTL', 'FNF'])]])
        rule_all_inputs.extend([[expand('output/qtl/{condition}_PEER_k{Nk}_batch_perm1Mb.txt', condition = ['CTL', 'FNF'], Nk = n) for n in range(10, 81, 10)]])
        pcQTL = ['output/qtl/{condition}_PC_batch_perm1Mb.txt']
        peerQTL = ['output/qtl/{condition}_PEER_k{Nk}_batch_perm1Mb.txt']
        rule_all_inputs.extend([[expand('output/qtl/{condition}_PC_batch_perm1Mb_FDR.txt', condition = ['CTL', 'FNF'])]])
        rule_all_inputs.extend([[expand('output/qtl/{condition}_PEER_k{Nk}_batch_perm1Mb_FDR.txt', condition = ['CTL', 'FNF'], Nk = n) for n in range(10, 81, 10)]])
        pcMultipleTesting = ['output/qtl/{condition}_PC_batch_perm1Mb_FDR.txt']
        peerMultipleTesting = ['output/qtl/{condition}_PEER_k{Nk}_batch_perm1Mb_FDR.txt']

        peer_eGene = [expand('output/qtl/{{condition}}_PEER_k{Nk}_batch_perm1Mb_FDR.txt', Nk = n) for n in range(10, 81, 10)]

        rule_all_inputs.extend([[expand('output/summary/{condition}_batch.txt', condition = ['CTL', 'FNF'])]])
        eGene = ['output/summary/{condition}_batch.txt']
    else:
        rule_all_inputs.extend([[expand('output/covar/{condition}_PC_covar.txt', condition = ['CTL', 'FNF'])]])
        rule_all_inputs.extend([[expand('output/covar/{condition}_PEER_k{Nk}_covar.txt', condition = ['CTL', 'FNF'], Nk = n) for n in range(10, 81, 10)]])
        pcCov = ['output/covar/{condition}_PC_covar.txt']
        peerCov = ['output/covar/{condition}_PEER_k{Nk}_covar.txt']
        rule_all_inputs.extend([[expand('output/qtl/{condition}_PC_perm1Mb.txt', condition = ['CTL', 'FNF'])]])
        rule_all_inputs.extend([[expand('output/qtl/{condition}_PEER_k{Nk}_perm1Mb.txt', condition = ['CTL', 'FNF'], Nk = n) for n in range(10, 81, 10)]])
        pcQTL = ['output/qtl/{condition}_PC_perm1Mb.txt']
        peerQTL = ['output/qtl/{condition}_PEER_k{Nk}_perm1Mb.txt']
        rule_all_inputs.extend([[expand('output/qtl/{condition}_PC_perm1Mb_FDR.txt', condition = ['CTL', 'FNF'])]])
        rule_all_inputs.extend([[expand('output/qtl/{condition}_PEER_k{Nk}_perm1Mb_FDR.txt', condition = ['CTL', 'FNF'], Nk = n) for n in range(10, 81, 10)]])
        pcMultipleTesting = ['output/qtl/{condition}_PC_perm1Mb_FDR.txt']
        peerMultipleTesting = ['output/qtl/{condition}_PEER_k{Nk}_perm1Mb_FDR.txt']

        peer_eGene = [expand('output/qtl/{{condition}}_PEER_k{Nk}_perm1Mb_FDR.txt', Nk = n) for n in range(10, 81, 10)]

        rule_all_inputs.extend([[expand('output/summary/{condition}_none.txt', condition = ['CTL', 'FNF'])]])
        eGene = ['output/summary/{condition}_none.txt']

## Define rules
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


# rule mbv:
#     input:
#         bam = rules.align.output,
#         index = rules.index.output,
#         vcf = rules.renameVCFdonors.output.vcf
#     output:
#         'output/mbv/{group}.bamstat.txt'
#     params:
#         version = config['QTLToolsVersion']
#     shell:
#         """
#         module load qtltools/{params.version}
#         QTLtools mbv --bam {input.bam} --vcf {input.vcf} --out output/mbv/{wildcards.group}.bamstat.txt
#         """

# # rule mbvParse:

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

rule getPCs:
    input:
        lambda wildcards: ['output/normquant/{condition}_CPMadjTMM_invNorm.bed.gz'.format(condition=wildcards.condition)]
    output:
        pcaCovar = 'output/covar/{condition}_PC.csv',
        pcaCor = 'output/pca/{condition}_PC.csv'
    params:
        version = config['Rversion'],
        samplesheet = config['samplesheet']
    log:
        out = 'output/logs/{condition}_getPCs.out',
        err = 'output/logs/{condition}_getPCs.err'
    shell:
        """
        module load r/{params.version}
        Rscript scripts/getPCs.R {params.samplesheet} {input} {wildcards.condition} 1> {log.out} 2> {log.err}
        """

rule getPCs_all:
    input:
        'output/normquant/ALL_CPMadjTMM_invNorm.bed.gz'
    output:
        'output/covar/ALL_PC.csv'
    params:
        version = config['Rversion'],
        samplesheet = config['samplesheet']
    log:
        out = 'output/logs/ALL_getPCs.out',
        err = 'output/logs/ALL_getPCs.err'
    shell:
        """
        module load r/{params.version}
        Rscript scripts/getPCs_all.R {params.samplesheet} {input} 1> {log.out} 2> {log.err}
        """

rule PCAcorrelation:
    input: 
        rules.getPCs.output.pcaCor
    output:
        'output/plots/{condition}_pca.pdf',
        'output/plots/{condition}_screeplot.pdf',
        'output/plots/{condition}_PC_correlation.png'
    params:
        version = config['Rversion'],
        samplesheet = config['samplesheet'],
        donorSamplesheet = config['donorSamplesheet']
    log:
        out = 'output/logs/{condition}_PCAcorrelation.out',
        err = 'output/logs/{condition}_PCAcorrelation.err'
    shell:
        """
        module load r/{params.version}
        Rscript scripts/PCAcorrelation.R {input} {params.donorSamplesheet} {params.samplesheet} {wildcards.condition} 1> {log.out} 2> {log.err}
        """

rule pcaKneedle:
    input:
        rules.getPCs.output.pcaCor
    output:
        'output/covar/{condition}_PCkneedle.txt'
    log:
        out = 'output/logs/{condition}_pcaKneedle.out',
        err = 'output/logs/{condition}_pcaKneedle.err'
    run:
        PC = pd.read_csv(str(input))

        # Grab sdev
        sdev = PC.loc[PC["Donor"] == 'sdev'].to_numpy()

        # Get rid of 'sdev' in first index
        sdev = np.delete(sdev, 0)

        # Calculate variance explained by each PC
        varExplained = np.square(sdev)/sum(np.square(sdev))

        # Generate vector for number of PCs
        x = range(1, len(PC.columns))

        # Calculate knee
        kneedle = KneeLocator(x, varExplained, curve = "convex", direction = "decreasing")

        # Write knee to file
        numPCs = kneedle.knee
        file = open(str(output), 'w')
        file.write(str(numPCs))
        file.close()


rule getPEER:
    input:
        lambda wildcards: ['output/normquant/{condition}_CPMadjTMM_invNorm.bed.gz'.format(condition=wildcards.condition)]
    output:
        [expand('output/covar/{{condition}}_PEERfactors_k{Nk}.txt', Nk = n) for n in range(10, 81, 10)]
    log:
        out = 'output/logs/{condition}_getPEER.out',
        err = 'output/logs/{condition}_getPEER.err'
    shell:
        """
        module load r/4.2.0
        Rscript scripts/PEERfactors.R {input} {wildcards.condition} 1> {log.out} 2> {log.err}
        """

rule checkPEER:
    input:
        [expand("output/covar/{{condition}}_PEERfactors_k{Nk}.txt", Nk = n) for n in range(10, 81, 10)]
    output:
        'output/covar/{condition}_validPEER.txt'
    log:
        out = 'output/logs/{condition}_checkPEER.out',
        err = 'output/logs/{condition}_checkPEER.err'
    params:
        condition = lambda wildcards: wildcards.condition
    run:
        for file in input:
            peerFactors = pd.read_csv(str(file))
            peerFactors = peerFactors.drop('Donor', axis = 1)
            if not peerFactors.isnull().values.all():
                index = input.index(file)
                Nk = range(10, 81, 10)
                # Write Nk as a valid PEER factor
                outFile = open(str(output), 'a')
                outFile.write(str(Nk[index])+"\n")
                outFile.close()

rule PEERcorrelation:
    input: 
        peer = lambda wildcards: 'output/covar/{condition}_PEERfactors_k{Nk}.txt'.format(condition=wildcards.condition, Nk=wildcards.Nk),
        check = rules.checkPEER.output
    output:
        'output/plots/{condition}_PEER{Nk}_correlation.png'
    params:
        version = config['Rversion'],
        samplesheet = config['samplesheet'],
        donorSamplesheet = config['donorSamplesheet']
    log:
        out = 'output/logs/{condition}_PEERk{Nk}_correlation.out',
        err = 'output/logs/{condition}_PEERk{Nk}_correlation.err'
    shell:
        """
        module load r/{params.version}
        Rscript scripts/PEERcorrelation.R {input.peer} {params.donorSamplesheet} {params.samplesheet} {wildcards.condition} {wildcards.Nk} {input.check} 1> {log.out} 2> {log.err}
        """

rule PC_eQTL:
    input:
        vcf = rules.renameVCFdonors.output.vcf,
        vcfIndex = rules.renameVCFdonors.output.index,
        bed = rules.indexQuant.output.bed,
        bedIndex = rules.indexQuant.output.index,
        cov = pcCov
    output:
        pcQTL
    params:
        version = config['QTLToolsVersion']
    log:
        out = 'output/logs/{condition}_PCeQTL.out',
        err = 'output/logs/{condition}_PCeQTL.err'
    shell:
        """
        module load qtltools/{params.version}
        QTLtools cis --vcf {input.vcf} --bed {input.bed} --cov {input.cov} --permute 1000 --window 1000000 --out {output} 1> {log.out} 2> {log.err}
        """

rule PEER_eQTL:
    input:
        vcf = rules.renameVCFdonors.output.vcf,
        vcfIndex = rules.renameVCFdonors.output.index,
        bed = rules.indexQuant.output.bed,
        bedIndex = rules.indexQuant.output.index,
        cov = peerCov,
        check = lambda wildcards: 'output/covar/{condition}_validPEER.txt'.format(condition=wildcards.condition)
    output:
        peerQTL
    params:
        version = config['QTLToolsVersion']
    log:
        out = 'output/logs/{condition}_PEERk{Nk}_eQTL.out',
        err = 'output/logs/{condition}_PEERk{Nk}_eQTL.err'
    shell:
        """
        module load qtltools/{params.version}
        if grep -Fx {wildcards.Nk} {input.check}
        then
            QTLtools cis --vcf {input.vcf} --bed {input.bed} --cov {input.cov} --permute 1000 --window 1000000 --out {output} 1> {log.out} 2> {log.err}
        else
            touch {output}
        fi
        """
    
rule PC_multipleTesting:
    input:
        qtlResult = pcQTL,
        geneInfo = rules.indexQuant.output.bed
    output:
        pcMultipleTesting
    params:
        version = config['Rversion']
    log:
        out = 'output/logs/{condition}_PCeQTL_multipleTesting.out',
        err = 'output/logs/{condition}_PCeQTL_multipleTesting.err'
    shell:
        """
        module load r/{params.version}
        Rscript scripts/correctQTLs.R {input.qtlResult} {input.geneInfo} {output} 1> {log.out} 2> {log.err}
        """

rule PEER_multipleTesting:
    input:
        qtlResult = peerQTL,
        geneInfo = rules.indexQuant.output.bed,
        check = lambda wildcards: 'output/covar/{condition}_validPEER.txt'.format(condition=wildcards.condition)
    output:
        peerMultipleTesting
    params:
        version = config['Rversion']
    log:
        out = 'output/logs/{condition}_PEERk{Nk}_eQTL_multipleTesting.out',
        err = 'output/logs/{condition}_PEERk{Nk}_eQTL_multipleTesting.err',
    shell:
        """
        module load r/{params.version}
        if grep -Fx {wildcards.Nk} {input.check}
        then
            Rscript scripts/correctQTLs.R {input.qtlResult} {input.geneInfo} {output} 1> {log.out} 2> {log.err}
        else
            touch {output}
        fi
        """

rule compareSig_eGenes:
    input:
        pcMultipleTesting,
        peer_eGene
    output:
        eGene
    params:
        version = config['Rversion'],
        FDRthreshold = config['FDRthreshold']
    log:
        out = 'output/logs/{condition}_compareSig_eGenes.out',
        err = 'output/logs/{condition}_compareSig_eGenes.err'
    shell:
        """
        module load r/{params.version}
        Rscript scripts/compareSig_eGenes.R {params.FDRthreshold} {output} {input} 1> {log.out} 2> {log.err}
        """