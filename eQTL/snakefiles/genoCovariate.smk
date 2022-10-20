if config['batchCovar'] == 'TRUE':
    pccovarOutput = ['output/covar/{condition}_PC_genoPC_batch_covar.txt']
    peercovarOutput = ['output/covar/{condition}_PEER_k{Nk}_genoPC_batch_covar.txt']
else:
    pccovarOutput = ['output/covar/{condition}_PC_genoPC_covar.txt'],
    peercovarOutput = ['output/covar/{condition}_PEER_k{Nk}_genoPC_covar.txt']

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
        batchCovar = config['batchCovar']
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

rule makePEERcovar_geno:
    input:
        peer = lambda wildcards: 'output/covar/{condition}_PEERfactors_k{Nk}.txt'.format(condition=wildcards.condition, Nk=wildcards.Nk),
        check = lambda wildcards: 'output/covar/{condition}_validPEER.txt'.format(condition=wildcards.condition),
        genoPC = rules.genoPCA.output.pcs,
        numgenoPCs = rules.genoPCkneedle.output
    output:
        peercovarOutput
    params:
        version = config['Rversion'],
        donorSamplesheet = config['donorSamplesheet'],
        samplesheet = config['samplesheet'],
        batchCovar = config['batchCovar']
    log:
        out = 'output/logs/{condition}_makePEERcovar_geno_k{Nk}.out',
        err = 'output/logs/{condition}_makePEERcovar_geno_k{Nk}.err'
    shell:
        """
        module load r/{params.version}
        numgenoPCs=`cat {input.numgenoPCs}`
        Rscript scripts/formatPEERcovariates_geno.R {input.peer} {params.donorSamplesheet} {params.samplesheet} {wildcards.condition} {wildcards.Nk} {input.check} {params.batchCovar} {input.genoPC} ${{numgenoPCs}} 1> {log.out} 2> {log.err}
        """