rule PEER_kneedle:
    input:
        [expand('output/covar/{{condition}}_PEERfactors_k{Nk}_variance.txt', Nk = Nk)]
    output:
        'output/covar/{condition}_PEERkneedle.txt'
    params:
        Nk = config['PEERfactors']
    log:
        out = 'output/logs/{condition}_PEERkneedle.out',
        err = 'output/logs/{condition}_PEERkneedle.err'
    run:
        var = pd.read_csv(str(input))["V1"].to_numpy()
        
        # Generate vector for number of PEER factors
        x = range(1, params.Nk + 1)

        # Calculate knee
        kneedle = KneeLocator(x, var, curve = "convex", direction = "decreasing")

        # Write knee to file
        numPEER = kneedle.knee
        file = open(str(output), 'w')
        file.write(str(numPEER))
        file.close()