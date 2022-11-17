batches = ['RNAKitBatch', 'RNASequencingBatch', 'genoBatch', 'DNAKitBatch']
if config['iteratePEER'] == 'TRUE':
    numPEER = 'None'
    peer = lambda wildcards: 'output/covar/{condition}_PEERfactors_k{Nk}.txt'.format(condition=wildcards.condition, Nk=wildcards.Nk)
    peercovarOutput = 'output/covar/{condition}_PEER_k{Nk}_genoPC'
    covar_logout = 'output/logs/{condition}_makePEERcovar_geno_k{Nk}.out'
    covar_logerr = 'output/logs/{condition}_makePEERcovar_geno_k{Nk}.err'
else:
    numPEER = lambda wildcards: 'output/covar/{condition}_PEERkneedle.txt'.format(condition=wildcards.condition)
    peer = lambda wildcards: 'output/covar/{condition}_PEERfactors_k{Nk}.txt'.format(condition=wildcards.condition, Nk=Nk)
    peercovarOutput = 'output/covar/{condition}_PEER_kneedle_genoPC'
    covar_logout = 'output/logs/{condition}_makePEERcovar_geno_kneedle.out'
    covar_logerr = 'output/logs/{condition}_makePEERcovar_geno_kneedle.err'

for b in batches:
    b_include = config[b]
    if b_include == "TRUE":
        peercovarOutput += '_{}'.format(b)

peercovarOutput += ".txt"

rule genoPCA:
    input:
        vcf = 'output/vcf/' + vcf_prefix + '_qtl.recode.vcf.gz',
        index = 'output/vcf/' + vcf_prefix + '_qtl.recode.vcf.gz.tbi'
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

rule makePEERcovar_geno:
    input:
        peer = peer,
        genoPC = rules.genoPCA.output.pcs,
        numgenoPCs = rules.genoPCkneedle.output
    output:
        peercovarOutput
    params:
        version = config['Rversion'],
        donorSamplesheet = config['donorSamplesheet'],
        dnaSamplesheet = config['dnaSamplesheet'],
        samplesheet = config['samplesheet'],
        RNAKitBatch = config['RNAKitBatch'],
        RNASequencingBatch = config['RNASequencingBatch'],
        genoBatch = config['genoBatch'],
        DNAKitBatch = config['DNAKitBatch'],
        numPEER = numPEER
    log:
        out = covar_logout,
        err = covar_logerr
    shell:
        """
        module load r/{params.version}
        numgenoPCs=`cat {input.numgenoPCs}`

        Rscript scripts/formatPEERcovariates_geno.R {input.peer} {params.numPEER} {params.donorSamplesheet} {params.dnaSamplesheet} {params.genoBatch} {params.DNAKitBatch} {params.samplesheet} {params.RNAKitBatch} {params.RNASequencingBatch} {wildcards.condition} {input.genoPC} ${{numgenoPCs}} {output} 1> {log.out} 2> {log.err}
        """