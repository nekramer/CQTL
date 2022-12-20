rule checkPEER:
    input:
        [expand("output/covar/{{condition}}_PEERfactors_k{Nk}.txt", Nk = Nk)]
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