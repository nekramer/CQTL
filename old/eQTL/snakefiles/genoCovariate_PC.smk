rule makePCcovar_geno:
    input:
        seqPC = lambda wildcards: ['output/covar/{condition}_PC.csv'.format(condition=wildcards.condition)],
        numPCs = lambda wildcards: ['output/covar/{condition}_PCkneedle.txt'.format(condition=wildcards.condition)],
        genoPC = rules.genoPCA.output.pcs,
        numgenoPCs = rules.genoPCkneedle.output
    output:
        pccovarOutput
    params:
        version = config['Rversion'],
        donorSamplesheet = config['donorSamplesheet'],
        samplesheet = config['samplesheet'],
        RNAKitBatch = config['RNAKitBatch']
    log:
        out = 'output/logs/{condition}_makePCcovar_geno.out',
        err = 'output/logs/{condition}_makePCcovar_geno.err'
    shell:
        """
        module load r/{params.version}
        numPCs=`cat {input.numPCs}`
        numgenoPCs=`cat {input.numgenoPCs}`
        Rscript scripts/formatPCcovariates_geno.R {input.seqPC} {params.donorSamplesheet} {params.samplesheet} ${{numPCs}} {input.genoPC} ${{numgenoPCs}} {params.batchCovar} {wildcards.condition} 1> {log.out} 2> {log.err}
        """