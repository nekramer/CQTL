rule chromContig:
    input:
        vcf
    output:
        temp('output/vcf/' + vcf_prefix + '_renameContig.vcf')
    threads: 1
    log:
        err = "output/vcf/logs/renameContig.err"
    params:
        chromNames = config['chromNames']
    shell:
        """
        module load samtools
        bcftools annotate --rename-chrs {params.chromNames} --threads {threads} -o {output} {input} 2> {log.err}
        """

rule renameChroms:
    input:
        rules.chromContig.output
    output:
        'output/vcf/' + vcf_prefix + '_renameChroms.vcf'
    threads: 1
    log:
        err = "output/vcf/logs/renameChroms.err"
    params:
        chromNames = config['chromNames']
    shell:
        """
        ## Look in second column to see if "chr" should be included or not
        if awk '{{print $2}}' {params.chromNames} | grep "chr"
        then
            # Add "chr" to ID column if necessary
            awk '{{if($0 !~ /^#/ && $3 !~ /^chr/) $3="chr"$3; print}}' {input} > {output} 2> {log.err}
        else
            # Remove "chr" from ID column
            awk '{{gsub(/\chr/, "")}}1' {input} > {output} 2> {log.err}
        fi
        """

rule zipVCF1:
    input:
        rules.renameChroms.output
    output:
        'output/vcf/' + vcf_prefix + '_renameChroms.vcf.gz'
    threads: 4
    log:
        err = "output/vcf/logs/zipVCF1.err"
    shell:
        """
        module load samtools
        bgzip --threads {threads} {input} 2> {log.err}
        """  

rule indexVCF1:
    input:
        rules.zipVCF1.output
    output:
        'output/vcf/' + vcf_prefix + '_renameChroms.vcf.gz.tbi'
    threads: 1
    log:
        out = "output/vcf/logs/indexVCF1.out",
        err = "output/vcf/logs/indexVCF1.err"
    shell:
        """
        module load gatk/4.1.7.0
        gatk IndexFeatureFile -I {input} 2> {log.err} 1> {log.out}
        """    

rule updateConfig:
    input:
        v = rules.zipVCF1.output,
        i = rules.indexVCF1.output
    output:
        v = 'output/vcf/' + vcf_prefix + '_newcontig.vcf.gz',
        i = 'output/vcf/' + vcf_prefix + '_newcontig.vcf.gz.tbi'
    threads: 1
    log:
        out = "output/vcf/logs/updateConfig.out",
        err = "output/vcf/logs/updateConfig.err"
    params:
        sequence = config['sequence']
    shell:
        """
        module load gatk/4.1.7.0
        gatk UpdateVCFSequenceDictionary -V {input.v} --source-dictionary {params.sequence} --output {output.v} --replace=true 2> {log.err} 1> {log.out}
        """

# rule selectVariants:
#     input:
#         v = rules.updateConfig.output.v,
#         i = rules.updateConfig.output.i
#     output:
#         v = temp('output/vcf/' + vcf_prefix + '_biallelic.vcf.gz'),
#         i = temp('output/vcf/' + vcf_prefix + '_biallelic.vcf.gz.tbi')
#     threads: 1
#     log:
#         out = "output/vcf/logs/selectVariants.out",
#         err = "output/vcf/logs/selectVariants.err"
#     params:
#         sequence = config['sequence']
#     shell:
#         """
#         module load gatk/4.1.7.0
#         gatk SelectVariants --variant {input.v} -R {params.sequence} --select-type-to-include SNP -O {output.v} --restrict-alleles-to BIALLELIC 2> {log.err} 1> {log.out}
#         """

# rule removeDuplicates:
#     input:
#         v = rules.selectVariants.output.v,
#         i = rules.selectVariants.output.i
#     output:
#         'output/vcf/' + vcf_prefix + '_nodups_biallelic.vcf'
#     threads: 4
#     log:
#         err = "output/vcf/logs/removeDuplicates.err"
#     shell:
#         """
#         module load samtools
#         bcftools norm -d any --threads {threads} -o {output} {input.v} 2> {log.err}
#         """

# rule zipVCF2:
#     input:
#         rules.removeDuplicates.output
#     output:
#         temp('output/vcf/' + vcf_prefix + '_nodups_biallelic.vcf.gz')
#     threads: 4
#     log:
#         err = "output/vcf/logs/zipVCF2.err"
#     shell:
#         """
#         module load samtools
#         bgzip --threads {threads} {input} 2> {log.err}
#         """  

# rule indexVCF2:
#     input:
#         rules.zipVCF2.output
#     output:
#         temp('output/vcf/' + vcf_prefix + '_nodups_biallelic.vcf.gz.tbi')
#     threads: 1
#     log:
#         out = "output/vcf/logs/indexVCF2.out",
#         err = "output/vcf/logs/indexVCF2.err"
#     shell:
#         """
#         module load samtools
#         tabix -p vcf {input} 2> {log.err} 1> {log.out}
#         """