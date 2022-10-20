if config['batchCovar'] == 'TRUE':
    pccovarOutput = ['output/covar/{condition}_PC_batch_covar.txt']
    peercovarOutput = ['output/covar/{condition}_PEER_k{Nk}_batch_covar.txt']
else:
    pccovarOutput = ['output/covar/{condition}_PC_covar.txt']
    peercovarOutput = ['output/covar/{condition}_PEER_k{Nk}_covar.txt']

rule makePCcovar:
    input:
        PC = lambda wildcards: ['output/covar/{condition}_PC.csv'.format(condition=wildcards.condition)],
        numPCs = lambda wildcards: ['output/covar/{condition}_PCkneedle.txt'.format(condition=wildcards.condition)]
    output:
        pccovarOutput
    params:
        version = config['Rversion'],
        donorSamplesheet = config['donorSamplesheet'],
        samplesheet = config['samplesheet'],
        batchCovar = config['batchCovar']
    log:
        out = 'output/logs/{condition}_makePCcovar.out',
        err = 'output/logs/{condition}_makePCcovar.err'
    shell:
        """
        module load r/{params.version}
        numPCs=`cat {input.numPCs}`
        Rscript scripts/formatPCcovariates.R {input.PC} {params.donorSamplesheet} {params.samplesheet} ${{numPCs}} {wildcards.condition} {params.batchCovar} 1> {log.out} 2> {log.err}
        """

rule makePEERcovar:
    input:
        peer = lambda wildcards: 'output/covar/{condition}_PEERfactors_k{Nk}.txt'.format(condition=wildcards.condition, Nk=wildcards.Nk),
        check = lambda wildcards: 'output/covar/{condition}_validPEER.txt'.format(condition=wildcards.condition)
    output:
        peercovarOutput
    params:
        version = config['Rversion'],
        donorSamplesheet = config['donorSamplesheet'],
        samplesheet = config['samplesheet'],
        batchCovar = config['batchCovar']
    log:
        out = 'output/logs/{condition}_makePEERcovar_k{Nk}.out',
        err = 'output/logs/{condition}_makePEERcovar_k{Nk}.err'
    shell:
        """
        module load r/{params.version}
        Rscript scripts/formatPEERcovariates.R {input.peer} {params.donorSamplesheet} {params.samplesheet} {wildcards.condition} {wildcards.Nk} {input.check} {params.batchCovar} 1> {log.out} 2> {log.err}
        """