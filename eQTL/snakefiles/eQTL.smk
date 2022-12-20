#!/usr/bin/env python3

# Parse minor allele filters
filterNum = config['minorAllele'].split(":")[0]
filterType = config['minorAllele'].split(":")[1]
if filterType == 'freq':
    filterFlag = '-q'

elif filterType == 'count':
    filterFlag = '-c'


if config['iteratePEER'] == 'TRUE':
    factors = [expand('output/covar/{{condition}}_PEERfactors_k{Nk}.txt', Nk = n) for n in range(5, Nk + 1, config['iterateBy'])],
    var = [expand('output/covar/{{condition}}_PEERfactors_k{Nk}_variance.txt', Nk = n) for n in range(5, Nk + 1, config['iterateBy'])]
else:
    include: "../rules/PEER_kneedle.smk"
    factors = [expand('output/covar/{{condition}}_PEERfactors_k{Nk}.txt', Nk = Nk)],
    var = [expand('output/covar/{{condition}}_PEERfactors_k{Nk}_variance.txt', Nk = Nk)]
    
include: "../../rules/VCFprocessing.smk"

include: "../../rules/RNAprocessing.smk"

rule subsetVCFdonors:
    input:
        vcf = rules.updateConfig.output.v,
        index = rules.updateConfig.output.i
    output:
        vcf = 'output/vcf/' + vcf_prefix + '_subset.vcf'
    params:
        donors = ",".join(samples['Donor'].unique().tolist()),
        samtoolsVersion = config['samtoolsVersion'],
        pythonVersion = config['pythonVersion'],
        prefix = vcf_prefix
    log:
        pythonOut = 'output/logs/subsetVCFdonors_py.out',
        pythonErr = 'output/logs/subsetVCFdonors_py.err',
        bcftoolsErr = 'output/logs/subsetVCFdonors_bcftools.err'
    shell:
        """
        module load samtools/{params.samtoolsVersion}
        module load python/{params.pythonVersion}
        bcftools query -l {input.vcf} > donors.txt
        python3 scripts/subsetVCFdonors.py donors.txt {params.donors} 1> {log.pythonOut} 2> {log.pythonErr}
        bcftools view -S subset.txt {input.vcf} > output/vcf/{params.prefix}_subset.vcf 2> {log.bcftoolsErr}
        rm donors.txt
        rm subset.txt
        """

rule renameVCFdonors:
    input:
        vcf = rules.subsetVCFdonors.output.vcf
    output:
        vcf = 'output/vcf/' + vcf_prefix + '_rename.vcf.gz',
        index = 'output/vcf/' + vcf_prefix + '_rename.vcf.gz.tbi'
    params:
        donors = ",".join(samples['Donor'].unique().tolist()),
        samtoolsVersion = config['samtoolsVersion'],
        pythonVersion = config['pythonVersion'],
        prefix = vcf_prefix
    log:
        pythonOut = 'output/logs/renameVCFdonors_py.out',
        pythonErr = 'output/logs/renameVCFdonors_py.err',
        bcftoolsErr = 'output/logs/renameVCFdonors_bcftools.err'
    shell:
        """
        module load samtools/{params.samtoolsVersion}
        module load python/{params.pythonVersion}
        bcftools query -l {input.vcf} > donors.txt
        python3 scripts/renameVCFdonors.py donors.txt {params.donors} 1> {log.pythonOut} 2> {log.pythonErr}
        bcftools reheader -s samples.txt -o output/vcf/{params.prefix}_rename.vcf {input.vcf} 2> {log.bcftoolsErr}
        bgzip output/vcf/{params.prefix}_rename.vcf && tabix -p vcf {output.vcf}
        rm donors.txt
        rm samples.txt
        """

## Rule for MAF filtering and het num filtering on VCF
rule filterVCFvariants:
    input:
        vcf = rules.renameVCFdonors.output.vcf,
        index = rules.renameVCFdonors.output.index
    output:
        vcf = 'output/vcf/' + vcf_prefix + '_qtl.recode.vcf.gz',
        index = 'output/vcf/' + vcf_prefix + '_qtl.recode.vcf.gz.tbi'
    params:
        samtoolsVersion = config['samtoolsVersion'],
        gatkVersion = config['gatkVersion'],
        vcftoolsVersion = config['vcftoolsVersion'],
        filterFlag = filterFlag,
        filterNum = filterNum,
        minHets = config['minHets'],
        prefix = vcf_prefix
    log:
        mafFilter_err = 'output/logs/mafFilter.err',
        hetFilter_out = 'output/logs/hetFilter.out',
        hetFilter_err = 'output/logs/hetFilter.err',
        vcfFilter_out = 'output/logs/vcfFilter.out',
        vcfFilter_err = 'output/logs/vcfFilter.err'
        
    shell:
        """
        module load samtools/{params.samtoolsVersion}
        module load gatk/{params.gatkVersion}
        module load vcftools/{params.vcftoolsVersion}
        bcftools view {params.filterFlag} {filterNum}:minor {input.vcf} > output/vcf/{params.prefix}_minorAlleleFilter.vcf 2> {log.mafFilter_err}
        gatk VariantFiltration -V output/vcf/{params.prefix}_minorAlleleFilter.vcf -O output/vcf/{params.prefix}_hetFilter.vcf.gz --filter-expression "vc.getHetCount() >= {params.minHets}" --filter-name "minhets" 1> {log.hetFilter_out} 2> {log.hetFilter_err}
        vcftools --gzvcf output/vcf/{params.prefix}_hetFilter.vcf.gz --keep-filtered minhets --recode --out output/vcf/{params.prefix}_qtl 1> {log.vcfFilter_out} 2> {log.vcfFilter_err}
        bgzip output/vcf/{params.prefix}_qtl.recode.vcf && tabix -p vcf {output.vcf}
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
        "output/normquant/{condition}_CPMadjTMM_invNorm.bed"
    params:
        version = config['Rversion'],
        samplesheet = config['samplesheet']
    log:
        out = 'output/logs/{condition}_quantNorm.out',
        err = 'output/logs/{condition}_quantNorm.err'
    shell:
        """
        module load r/{params.version}
        Rscript scripts/quantNorm.R {params.samplesheet} {wildcards.condition} {input} 1> {log.out} 2> {log.err}
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
        factors = factors,
        var = var
    params:
        Nk = config['PEERfactors'],
        iteratePEER = config['iteratePEER'],
        iterateBy = config['iterateBy']
    log:
        out = 'output/logs/{condition}_getPEER.out',
        err = 'output/logs/{condition}_getPEER.err'
    shell:
        """
        module load r/4.2.0
        Rscript scripts/PEERfactors.R {input} {wildcards.condition} {params.Nk} {params.iteratePEER} {params.iterateBy} 1> {log.out} 2> {log.err}
        """