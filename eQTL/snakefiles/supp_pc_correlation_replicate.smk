#!/usr/bin/env python3
import pandas as pd

## Load config file
configfile: "config/config_pc_correlation.yaml"

## Read in samplesheet
samples = pd.read_csv(config["samplesheet"], sep = ",")
donors = pd.read_csv(config["donorSamplesheet"])
dna = pd.read_csv(config["dnaSamplesheet"])

## Get possible covariates from samplesheets
sample_covariates = list(samples.columns)
donor_covariates = list(donors.columns)
dna_covariates = list(dna.columns)
covariates = sample_covariates + donor_covariates + dna_covariates

## Convert samplesheet columns to strings
samples = samples.astype(str)

## Concatenate Sequencing_Directory to Read1 and Read2 for full read paths
samples['Read1'] = samples[['Sequencing_Directory', 'Read1']].apply(lambda row: os.path.join(*row), axis=1)
samples['Read2'] = samples[['Sequencing_Directory', 'Read2']].apply(lambda row: os.path.join(*row), axis=1)

## Group Seq_Reps
samples['id'] = samples[['Proj', 'Donor']].agg('_'.join, axis=1) + '_R_' + samples[['Condition', 'Time', 'Tech_Rep']].agg('_'.join, axis=1)

## Extract grouped read1 and read2
read1 = samples.groupby(['id'])['Read1'].apply(list).to_dict()

rule all:
    input:
        [expand("output/mbv/{group}.bamstat.txt", group = key) for key in read1], 
        [expand('output/mbv/{group}_concordance.txt', group = key) for key in read1],
        'output/mbv/mbvSummary.txt',
        [expand("output/pca/{condition}_PCs.csv", condition = ['ALL', 'CTL', 'FNF'])],
        [expand("output/pca/{condition}_PC_varExplained.csv", condition = ['ALL', 'CTL', 'FNF'])],
        [expand("output/pca/{condition}_PCs_RNAextractionKitBatch_corrected.csv", condition = ['ALL', 'CTL', 'FNF'])],
        [expand("output/pca/{condition}_PC_RNAextractionKitBatch_corrected_varExplained.csv", condition = ['ALL', 'CTL', 'FNF'])],
        [expand("output/pca/{condition}_PCs_SequencingBatch_corrected.csv", condition = ['ALL', 'CTL', 'FNF'])],
        [expand("output/pca/{condition}_PC_SequencingBatch_corrected_varExplained.csv", condition = ['ALL', 'CTL', 'FNF'])],
        [expand("output/pca/{condition}_PCs_DNAReagentBatch_corrected.csv", condition = ['ALL', 'CTL', 'FNF'])],
        [expand("output/pca/{condition}_PC_DNAReagentBatch_corrected_varExplained.csv", condition = ['ALL', 'CTL', 'FNF'])],
        [expand("output/pca/{condition}_PCs_GenotypingBatch_corrected.csv", condition = ['ALL', 'CTL', 'FNF'])],
        [expand("output/pca/{condition}_PC_GenotypingBatch_corrected_varExplained.csv", condition = ['ALL', 'CTL', 'FNF'])],
        [expand("output/pca/{condition}_PCs_FragmentBatch_corrected.csv", condition = ['ALL', 'CTL', 'FNF'])],
        [expand("output/pca/{condition}_PC_FragmentBatch_corrected_varExplained.csv", condition = ['ALL', 'CTL', 'FNF'])],
        [expand("output/pca/{condition}_correlations.csv", condition = ['ALL', 'CTL', 'FNF'])],
        [expand("output/pca/{condition}_correlations_pval.csv", condition = ['ALL', 'CTL', 'FNF'])],
        [expand("output/peer/{condition}_PEERcorrelations.csv", condition = ['ALL', 'CTL', 'FNF'])],
        [expand("output/peer/{condition}_PEERcorrelations_pval.csv", condition = ['ALL', 'CTL', 'FNF'])],
        [expand('output/pca/{condition}_covarPC.RData', condition = ['ALL', 'CTL', 'FNF'])],
        [expand('output/pca/{condition}_screeplot.pdf', condition = ['ALL', 'CTL', 'FNF'])],
        [expand('output/pca/{condition}_correlation.pdf', condition = ['ALL', 'CTL', 'FNF'])],
        [expand('output/peer/{condition}_PEERcorrelation.pdf', condition = ['ALL', 'CTL', 'FNF'])],
        [expand('output/GTEx/{condition}_eGene_GTEx_pi1.csv', condition = ['CTL', 'FNF'])],
        [expand('output/GTEx/{condition}_eGene_GTEx_downsampled_percOverlap.csv', condition = ['CTL', 'FNF'])],
        'output/GTEx/GTEx_eGene_pi1.pdf',
        'output/GTEx/GTEx_eGene_percentOverlap.pdf',
        [expand('output/pca/{condition}_RNAextractionKitBatch_corrected_PCplots.RData', condition = ['ALL', 'CTL', 'FNF'])],
        [expand('output/pca/{condition}_SequencingBatch_corrected_PCplots.RData', condition = ['ALL', 'CTL', 'FNF'])],
        [expand('output/pca/{condition}_DNAReagentBatch_corrected_PCplots.RData', condition = ['ALL', 'CTL', 'FNF'])],
        [expand('output/pca/{condition}_GenotypingBatch_corrected_PCplots.RData', condition = ['ALL', 'CTL', 'FNF'])],
        [expand('output/pca/{condition}_FragmentBatch_corrected_PCplots.RData', condition = ['ALL', 'CTL', 'FNF'])],
        [expand('output/pca/{condition}_PCs_RNAextractionKit_DNAReagent_corrected.csv', condition = ['ALL', 'CTL', 'FNF'])],
        [expand('output/pca/{condition}_PC_RNAextractionKit_DNAReagent_corrected_varExplained.csv', condition = ['ALL', 'CTL', 'FNF'])],
        [expand('output/pca/{condition}_RNAextractionKit_DNAReagent_corrected_correlations.csv', condition = ['ALL', 'CTL', 'FNF'])],
        [expand('output/pca/{condition}_RNAextractionKit_DNAReagent_corrected_correlations_pval.csv', condition = ['ALL', 'CTL', 'FNF'])],
        [expand('output/pca/{condition}_RNAextractionKitBatch_corrected_correlations.csv', condition = ['ALL', 'CTL', 'FNF'])],
        [expand('output/pca/{condition}_RNAextractionKitBatch_corrected_correlations_pval.csv', condition = ['ALL', 'CTL', 'FNF'])],
        [expand('output/pca/{condition}_SequencingBatch_corrected_correlations.csv', condition = ['ALL', 'CTL', 'FNF'])],
        [expand('output/pca/{condition}_SequencingBatch_corrected_correlations_pval.csv', condition = ['ALL', 'CTL', 'FNF'])],
        [expand('output/pca/{condition}_DNAReagentBatch_corrected_correlations.csv', condition = ['ALL', 'CTL', 'FNF'])],
        [expand('output/pca/{condition}_DNAReagentBatch_corrected_correlations_pval.csv', condition = ['ALL', 'CTL', 'FNF'])],
        [expand('output/pca/{condition}_GenotypingBatch_corrected_correlations.csv', condition = ['ALL', 'CTL', 'FNF'])],
        [expand('output/pca/{condition}_GenotypingBatch_corrected_correlations_pval.csv', condition = ['ALL', 'CTL', 'FNF'])],
        [expand('output/pca/{condition}_FragmentBatch_corrected_correlations.csv', condition = ['ALL', 'CTL', 'FNF'])],
        [expand('output/pca/{condition}_FragmentBatch_corrected_correlations_pval.csv', condition = ['ALL', 'CTL', 'FNF'])],
        [expand('output/pca/{condition}_RNAextractionKitBatch_corrected_correlation.pdf', condition = ['ALL', 'CTL', 'FNF'])],
        [expand('output/pca/{condition}_SequencingBatch_corrected_correlation.pdf', condition = ['ALL', 'CTL', 'FNF'])],
        [expand('output/pca/{condition}_DNAReagentBatch_corrected_correlation.pdf', condition = ['ALL', 'CTL', 'FNF'])],
        [expand('output/pca/{condition}_GenotypingBatch_corrected_correlation.pdf', condition = ['ALL', 'CTL', 'FNF'])],
        [expand('output/pca/{condition}_FragmentBatch_corrected_correlation.pdf', condition = ['ALL', 'CTL', 'FNF'])],
        [expand('output/pca/{condition}_RNAextractionKit_DNAReagent_corrected_PCplots.RData', condition = ['ALL', 'CTL', 'FNF'])],
        [expand('output/pca/{condition}_RNAextractionKit_DNAReagent_corrected_correlation.pdf', condition = ['ALL', 'CTL', 'FNF'])],
        'output/peer/PEER_eGene_comparisons.pdf'
        
rule mbv:
    input:
        bam = lambda wildcards: [config['data'] + '/{group}/align/{group}.Aligned.sortedByCoord.out.bam'.format(group = wildcards.group)],
        index = lambda wildcards: [config['data'] + '/{group}/align/{group}.Aligned.sortedByCoord.out.bam.bai'.format(group = wildcards.group)],
        vcf = config['vcf']
    output:
        'output/mbv/{group}.bamstat.txt'
    params:
        version = config['QTLToolsVersion']
    shell:
        """
        module load qtltools/{params.version}
        QTLtools mbv --bam {input.bam} --vcf {input.vcf} --out output/mbv/{wildcards.group}.bamstat.txt
        """

rule mbvParse:
    input:
        rules.mbv.output
    output:
        'output/mbv/{group}_concordance.txt'
    params:
        version = config['Rversion']
    log:
        out = 'output/logs/{group}_mbvParse.out',
        err = 'output/logs/{group}_mbvParse.err'
    shell:
        """
        module load r/{params.version}
        Rscript scripts/supp_pc_correlation_comparison/mbv_parsing.R {input} {wildcards.group} {output} 1> {log.out} 2> {log.err}
        """

rule combine_mbvParse:
    input:
        [expand('output/mbv/{group}_concordance.txt', group = key) for key in read1]
    output:
        'output/mbv/mbvSummary.txt'
    shell:
        """
        head -n 1 {input[1]} > output/mbv/mbvSummary.txt
        tail -n +2 -q {input} >> output/mbv/mbvSummary.txt
        """

rule normalizeIndels:
    input:
        config['vcf']
    output:
        vcf = 'output/supp/normalized_indels.vcf.gz',
        vcfi = 'output/supp/normalized_indelx.vcf.gz.tbi'
    params:
        version = config['samtoolsVersion']
    shell:
        """
        module load samtools/{params.version}
        bcftools norm -d all {input} -o output/supp/normalized_indels.vcf.gz
        tabix -p vcf output/supp/normalized_indels.vcf.gz
        """

rule generate_verifybamid_resourceFiles:
    input:
        vcf = rules.normalizeIndels.output.vcf,
        vcfi = rules.normalizeIndels.output.vcfi
    output:
        UD = 'output/supp/normalized_indels.vcf.gz.UD',
        mu = 'output/supp/normalized_indels.vcf.gz.mu',
        bed = 'output/supp/normalized_indels.vcf.gz.bed'
    params:
        version = config['pythonVersion'],
        sequence = config['sequence']
    log:
        out = 'output/logs/generate_verifybamid_resourceFiles.out',
        err = 'output/logs/generate_verifybamid_resourceFiles.err'
    shell:
        """
        module load python/{params.version}
        verifybamid2 --RefVCF {input.vcf} --Reference {params.sequence} 1> {log.out} 2> {log.err}
        """

rule verifybamid:
    input:
        bam = lambda wildcards: [config['data'] + '/{group}/align/{group}.Aligned.sortedByCoord.out.bam'.format(group = wildcards.group)],
        index = lambda wildcards: [config['data'] + '/{group}/align/{group}.Aligned.sortedByCoord.out.bam.bai'.format(group = wildcards.group)],
        UD = rules.generate_verifybamid_resourceFiles.output.UD,
        mu = rules.generate_verifybamid_resourceFiles.output.mu,
        bed = rules.generate_verifybamid_resourceFiles.output.bed
    output:

    params:
        version = config['pythonVersion'],
        sequence = config['sequence']
    log:
        out = 'output/logs/verifybamid.out',
        err = 'output/logs/verifybamid.err'
    shell:
        """
        module load python/{params.version}
        verifybamid2 --BamFile {input.bam} --SVDPrefix output/supp/normalized_indels.vcf.gz --Reference {params.sequence} --DisableSanityCheck 1> {log.out} 2> {log.err}
        """

rule getPCs:
    input:
        #lambda wildcards: [config['data'] + '/normquant/{condition}_CPMadjTMM_invNorm.bed.gz'.format(condition = wildcards.condition)]
        config['data'] + '/normquant/{condition}_CPMadjTMM_invNorm.bed.gz'
    output:
        pcs = 'output/pca/{condition}_PCs.csv',
        varExplained = 'output/pca/{condition}_PC_varExplained.csv'
    params:
        samplesheet = config['samplesheet'],
        dnaSamplesheet = config['dnaSamplesheet'],
        donorSamplesheet = config['donorSamplesheet'],
        version = config['Rversion']
    log:
        out = 'output/logs/{condition}_getPCs.out',
        err = 'output/logs/{condition}_getPCs.err'
    shell:
        """
        module load r/{params.version}
        Rscript scripts/supp_pc_correlation_comparison/getPCs.R {input} {wildcards.condition} {params.samplesheet} {params.dnaSamplesheet} {params.donorSamplesheet} 1> {log.out} 2> {log.err}
        """

rule PCAcorrelation:
    input: 
        lambda wildcards: ['output/pca/{condition}_PCs.csv'.format(condition = wildcards.condition)]
    output:
        pcCors = 'output/pca/{condition}_correlations.csv',
        pcCor_pval = 'output/pca/{condition}_correlations_pval.csv'
    params:
        version = config['Rversion']
    log:
        out = 'output/logs/{condition}_PCAcorrelation.out',
        err = 'output/logs/{condition}_PCAcorrelation.err'
    shell:
        """
        module load r/{params.version}
        Rscript scripts/supp_pc_correlation_comparison/PCAcorrelation.R {input} {wildcards.condition} {output.pcCors} {output.pcCor_pval} 1> {log.out} 2> {log.err}
        """

rule PEERcorrelation:
    input: 
        lambda wildcards: [config['peer'] + '/{condition}_PEERfactors_k'.format(condition = wildcards.condition) + str(config['Nk']) + '.txt']
    output:
        peer_cors = 'output/peer/{condition}_PEERcorrelations.csv',
        peer_cors_pval = 'output/peer/{condition}_PEERcorrelations_pval.csv'
    params:
        version = config['Rversion'],
        samplesheet = config['samplesheet'],
        dnaSamplesheet = config['dnaSamplesheet'],
        donorSamplesheet = config['donorSamplesheet']
    log:
        out = 'output/logs/{condition}_PEERcorrelation.out',
        err = 'output/logs/{condition}_PEERcorrelation.err'
    shell:
        """
        module load r/{params.version}
        Rscript scripts/supp_pc_correlation_comparison/PEERcorrelation.R {input} {wildcards.condition} {params.samplesheet} {params.dnaSamplesheet} {params.donorSamplesheet} 1> {log.out} 2> {log.err}
        """

rule make_pc_corr_Plots:
    input:
        pcs = rules.getPCs.output.pcs,
        varExplained = rules.getPCs.output.varExplained,
        pcCors = rules.PCAcorrelation.output.pcCors,
        pcCor_pval = rules.PCAcorrelation.output.pcCor_pval,
        peer_cors = rules.PEERcorrelation.output.peer_cors,
        peer_cors_pval = rules.PEERcorrelation.output.peer_cors_pval
    output:
        'output/pca/{condition}_covarPC.RData',
        'output/pca/{condition}_screeplot.pdf',
        'output/pca/{condition}_correlation.pdf',
        'output/peer/{condition}_PEERcorrelation.pdf'
    params:
        version = config['Rversion']
    log:
        out = 'output/logs/{condition}_make_pc_corr_Plots.out',
        err = 'output/logs/{condition}_make_pc_corr_Plots.err'
    shell:
        """
        module load r/{params.version}
        Rscript scripts/supp_pc_correlation_comparison/pc_correlation_plotting.R {input.pcs} {input.varExplained} {input.pcCors} {input.pcCor_pval} {input.peer_cors} {input.peer_cors_pval} {wildcards.condition} 1> {log.out} 2> {log.err}
        """

rule GTEx_comparisons:
    input:
        lambda wildcards: [config['eGenes'] + '/{condition}_PEER_k'.format(condition = wildcards.condition) + str(config['Nk']) + '_genoPC_RNAKitBatch_perm1Mb_FDR.txt']
    output:
        pi1 = 'output/GTEx/{condition}_eGene_GTEx_pi1.csv',
        percOverlap = 'output/GTEx/{condition}_eGene_GTEx_downsampled_percOverlap.csv'
    params:
        version = config['Rversion'],
        gtex_egene = config['gtex_egene'],
        gtex_signif = config['gtex_signif'],
        FDRthreshold = config['FDRthreshold']
    log:
        out = 'output/logs/{condition}_GTEx_comparisons.out',
        err = 'output/logs/{condition}_GTEx_comparisons.err'   
    shell:
        """
        module load r/{params.version}
        Rscript scripts/supp_pc_correlation_comparison/eGene_comparisons_GTEx.R {input} {params.gtex_egene} {params.gtex_signif} {params.FDRthreshold} {wildcards.condition} 1> {log.out} 2> {log.err}
        """

rule plot_GTEx_comparisons:
    input:
        CTL_pi1 = 'output/GTEx/CTL_eGene_GTEx_pi1.csv',
        FNF_pi1 = 'output/GTEx/FNF_eGene_GTEx_pi1.csv',
        CTL_percOverlap = 'output/GTEx/CTL_eGene_GTEx_downsampled_percOverlap.csv',
        FNF_percOverlap = 'output/GTEx/FNF_eGene_GTEx_downsampled_percOverlap.csv'
    output:
        'output/GTEx/GTEx_eGene_pi1.pdf',
        'output/GTEx/GTEx_eGene_percentOverlap.pdf'
    params:
        version = config['Rversion']
    log:
        out = 'output/logs/plot_GTEx_comparisons.out',
        err = 'output/logs/plot_GTEx_comparisons.err'
    shell:
        """
        module load r/{params.version}
        Rscript scripts/supp_pc_correlation_comparison/plot_GTEx_comparisons.R {input.CTL_pi1} {input.FNF_pi1} {input.CTL_percOverlap} {input.FNF_percOverlap} 1> {log.out} 2> {log.err}
        """

rule getPCs_batchCorrection:
    input:
        lambda wildcards: [config['data'] + '/normquant/{condition}_CPMadjTMM_invNorm.bed.gz'.format(condition = wildcards.condition)]
    output:
        'output/pca/{condition}_PCs_RNAextractionKitBatch_corrected.csv',
        'output/pca/{condition}_PC_RNAextractionKitBatch_corrected_varExplained.csv',
        'output/pca/{condition}_PCs_SequencingBatch_corrected.csv',
        'output/pca/{condition}_PC_SequencingBatch_corrected_varExplained.csv',
        'output/pca/{condition}_PCs_DNAReagentBatch_corrected.csv',
        'output/pca/{condition}_PC_DNAReagentBatch_corrected_varExplained.csv',
        'output/pca/{condition}_PCs_GenotypingBatch_corrected.csv',
        'output/pca/{condition}_PC_GenotypingBatch_corrected_varExplained.csv',
        'output/pca/{condition}_PCs_FragmentBatch_corrected.csv',
        'output/pca/{condition}_PC_FragmentBatch_corrected_varExplained.csv'
    params:
        samplesheet = config['samplesheet'],
        dnaSamplesheet = config['dnaSamplesheet'],
        donorSamplesheet = config['donorSamplesheet'],
        version = config['Rversion']
    log:
        out = 'output/logs/{condition}_getPCs_batchCorrection.out',
        err = 'output/logs/{condition}_getPCs_batchCorrection.err'
    shell:
        """
        module load r/{params.version}
        for batch in RNAextractionKitBatch SequencingBatch DNAReagentBatch GenotypingBatch FragmentBatch
        do
            Rscript scripts/supp_pc_correlation_comparison/getPCs_batchCorrection.R {input} {wildcards.condition} ${{batch}} {params.samplesheet} {params.dnaSamplesheet} {params.donorSamplesheet} 1> {log.out} 2> {log.err}
        done
        """

rule PCA_batchCorrection_correlation:
    input: 
        'output/pca/{condition}_PCs_RNAextractionKitBatch_corrected.csv',
        'output/pca/{condition}_PCs_SequencingBatch_corrected.csv',
        'output/pca/{condition}_PCs_DNAReagentBatch_corrected.csv',
        'output/pca/{condition}_PCs_GenotypingBatch_corrected.csv',
        'output/pca/{condition}_PCs_FragmentBatch_corrected.csv'
    output:
        'output/pca/{condition}_RNAextractionKitBatch_corrected_correlations.csv',
        'output/pca/{condition}_RNAextractionKitBatch_corrected_correlations_pval.csv',
        'output/pca/{condition}_SequencingBatch_corrected_correlations.csv',
        'output/pca/{condition}_SequencingBatch_corrected_correlations_pval.csv',
        'output/pca/{condition}_DNAReagentBatch_corrected_correlations.csv',
        'output/pca/{condition}_DNAReagentBatch_corrected_correlations_pval.csv',
        'output/pca/{condition}_GenotypingBatch_corrected_correlations.csv',
        'output/pca/{condition}_GenotypingBatch_corrected_correlations_pval.csv',
        'output/pca/{condition}_FragmentBatch_corrected_correlations.csv',
        'output/pca/{condition}_FragmentBatch_corrected_correlations_pval.csv'
    params:
        version = config['Rversion']
    log:
        out = 'output/logs/{condition}_PCA_batchCorrection_correlation.out',
        err = 'output/logs/{condition}_PCA_batchCorrection_correlation.err'
    shell:
        """
        module load r/{params.version}
        for batch in RNAextractionKitBatch SequencingBatch DNAReagentBatch GenotypingBatch FragmentBatch
        do
            Rscript scripts/supp_pc_correlation_comparison/PCAcorrelation.R output/pca/{wildcards.condition}_PCs_${{batch}}_corrected.csv {wildcards.condition} output/pca/{wildcards.condition}_${{batch}}_corrected_correlations.csv output/pca/{wildcards.condition}_${{batch}}_corrected_correlations_pval.csv 1> {log.out} 2> {log.err}
        done
        """

rule make_pc_corr_batchCorrection_plots:
    input:
        rules.getPCs_batchCorrection.output,
        rules.PCA_batchCorrection_correlation.output
    output:
        'output/pca/{condition}_RNAextractionKitBatch_corrected_PCplots.RData',
        'output/pca/{condition}_SequencingBatch_corrected_PCplots.RData',
        'output/pca/{condition}_DNAReagentBatch_corrected_PCplots.RData',
        'output/pca/{condition}_GenotypingBatch_corrected_PCplots.RData',
        'output/pca/{condition}_FragmentBatch_corrected_PCplots.RData',
        'output/pca/{condition}_RNAextractionKitBatch_corrected_correlation.pdf',
        'output/pca/{condition}_SequencingBatch_corrected_correlation.pdf',
        'output/pca/{condition}_DNAReagentBatch_corrected_correlation.pdf',
        'output/pca/{condition}_GenotypingBatch_corrected_correlation.pdf',
        'output/pca/{condition}_FragmentBatch_corrected_correlation.pdf'
    params:
        version = config['Rversion']
    log:
        out = 'output/logs/{condition}_plotPCs_batchCorrection.out',
        err = 'output/logs/{condition}_plotPCs_batchCorrection.err'
    shell:
        """
        module load r/{params.version}
        for batch in RNAextractionKitBatch SequencingBatch DNAReagentBatch GenotypingBatch FragmentBatch
        do
            Rscript scripts/supp_pc_correlation_comparison/pc_batchCorrection_plotting.R output/pca/{wildcards.condition}_PCs_${{batch}}_corrected.csv output/pca/{wildcards.condition}_PC_${{batch}}_corrected_varExplained.csv ${{batch}} output/pca/{wildcards.condition}_${{batch}}_corrected_correlations.csv output/pca/{wildcards.condition}_${{batch}}_corrected_correlations_pval.csv {wildcards.condition} 1> {log.out} 2> {log.err}
        done
        """

rule getPCs_multipleBatchCorrection:
    input:
        lambda wildcards: [config['data'] + '/normquant/{condition}_CPMadjTMM_invNorm.bed.gz'.format(condition = wildcards.condition)]
    output:
        'output/pca/{condition}_PCs_RNAextractionKit_DNAReagent_corrected.csv',
        'output/pca/{condition}_PC_RNAextractionKit_DNAReagent_corrected_varExplained.csv'
    params:
        version = config['Rversion'],
        samplesheet = config['samplesheet'],
        dnaSamplesheet = config['dnaSamplesheet'],
        donorSamplesheet = config['donorSamplesheet']
    log:
        out = 'output/logs/{condition}_getPCs_multipleBatchCorrection.out',
        err = 'output/logs/{condition}_getPCs_multipleBatchCorrection.err'
    shell:
        """
        module load r/{params.version}
        Rscript scripts/supp_pc_correlation_comparison/getPCs_multipleBatchCorrection.R {input} {wildcards.condition} {params.samplesheet} {params.dnaSamplesheet} {params.donorSamplesheet} 1> {log.out} 2> {log.err}
        """

rule PCA_multipleBatchCorrection_correlation:
    input: 
        'output/pca/{condition}_PCs_RNAextractionKit_DNAReagent_corrected.csv'
    output:
        pcCors = 'output/pca/{condition}_RNAextractionKit_DNAReagent_corrected_correlations.csv',
        pcCor_pval = 'output/pca/{condition}_RNAextractionKit_DNAReagent_corrected_correlations_pval.csv'
    params:
        version = config['Rversion']
    log:
        out = 'output/logs/{condition}_PCA_multipleBatchCorrection_correlation.out',
        err = 'output/logs/{condition}_PCA_multipleBatchCorrection_correlation.err'
    shell:
        """
        module load r/{params.version}
        Rscript scripts/supp_pc_correlation_comparison/PCAcorrelation.R {input} {wildcards.condition} {output.pcCors} {output.pcCor_pval} 1> {log.out} 2> {log.err}
        """

rule make_pc_corr_multipleBatchCorrection_plots:
    input:
        rules.getPCs_multipleBatchCorrection.output,
        rules.PCA_multipleBatchCorrection_correlation.output
    output:
        'output/pca/{condition}_RNAextractionKit_DNAReagent_corrected_PCplots.RData',
        'output/pca/{condition}_RNAextractionKit_DNAReagent_corrected_correlation.pdf'
    params:
        version = config['Rversion']
    log:
        out = 'output/logs/{condition}_make_pc_corr_multipleBatchCorrection.out',
        err = 'output/logs/{condition}_make_pc_corr_multipleBatchCorrection.err'
    shell:
        """
        module load r/{params.version}
        Rscript scripts/supp_pc_correlation_comparison/pc_batchCorrection_plotting.R output/pca/{wildcards.condition}_PCs_RNAextractionKit_DNAReagent_corrected.csv output/pca/{wildcards.condition}_PC_RNAextractionKit_DNAReagent_corrected_varExplained.csv RNAextractionKit_DNAReagent output/pca/{wildcards.condition}_RNAextractionKit_DNAReagent_corrected_correlations.csv output/pca/{wildcards.condition}_RNAextractionKit_DNAReagent_corrected_correlations_pval.csv {wildcards.condition} 1> {log.out} 2> {log.err}
        """

rule compare_PEER_eGenes:
    input:
        config['eGenes']
    output:
        'output/peer/PEER_eGene_comparisons.pdf'
    params:
        version = config['Rversion'],
        FDRthreshold = config['FDRthreshold']
    log:
        out = 'output/logs/compare_PEER_eGenes.out',
        err = 'output/logs/compare_PEER_eGenes.err'
    shell:
        """
        module load r/{params.version}
        Rscript scripts/supp_pc_correlation_comparison/plot_PEER_comparison.R {input} {params.FDRthreshold} 1> {log.out} 2> {log.err}
        """