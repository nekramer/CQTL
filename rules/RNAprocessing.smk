rule catR1:
    input:
        lambda wildcards: read1.get(wildcards.group)
    output:
        "output/{group}/fastq/{group}_R1.fastq.gz"
    threads: 1
    log:
        err = "output/{group}/logs/{group}_catR1.err"
    shell:
        """
        mkdir -p output/{wildcards.group}/fastq
        cat {input} > {output} 2> {log.err}
        """

rule catR2:
    input:
        lambda wildcards: read2.get(wildcards.group)
    output:
        "output/{group}/fastq/{group}_R2.fastq.gz"
    threads: 1
    log:
        err = "output/{group}/logs/{group}_catR2.err"
    shell:
        """
        mkdir -p output/{wildcards.group}/fastq
        cat {input} > {output} 2> {log.err}
        """

rule qc:
    input:
        R1 = lambda wildcards: ['output/{group}/fastq/{group}_R1.fastq.gz'.format(group=wildcards.group)],
        R2 = lambda wildcards: ['output/{group}/fastq/{group}_R2.fastq.gz'.format(group=wildcards.group)]
    output:
        zips = expand('output/qc/{{group}}_{R}_fastqc.zip', R=['R1', 'R2']),
        html = expand('output/qc/{{group}}_{R}_fastqc.html',R=['R1', 'R2'])
    threads: 2
    params:
        version = config['fastqcVersion']
    log:
        err = "output/{group}/logs/{group}_qc.err"
    shell:
        """
        module load fastqc/{params.version}
        mkdir -p output/qc
        fastqc -t {threads} -o output/qc/{wildcards.group} {input.R1} {input.R2} 2> {log.err}
        """

rule trim:
    input:
        R1 = lambda wildcards: ['output/{group}/fastq/{group}_R1.fastq.gz'.format(group=wildcards.group)],
        R2 = lambda wildcards: ['output/{group}/fastq/{group}_R2.fastq.gz'.format(group=wildcards.group)]
    output:
        trim1 = temp("output/{group}/trim/{group}_R1_val_1.fq.gz"),
        trim2 = temp("output/{group}/trim/{group}_R2_val_2.fq.gz"),
        report1 = temp("output/{group}/trim/{group}_R1_trimming_report.txt"),
        report2 = temp("output/{group}/trim/{group}_R2_trimming_report.txt")
    threads: 4
    params:
        version = config['trimgaloreVersion']
    log:
        err = "output/{group}/logs/{group}_trim.err"
    shell:
        """
        module load trim_galore/{params.version}
        module load python/3.6.6
        module load pigz
        mkdir -p output/{wildcards.group}/trim
        trim_galore -o output/{wildcards.group}/trim --cores {threads} --path_to_cutadapt /nas/longleaf/apps/cutadapt/2.9/venv/bin/cutadapt --paired {input.R1} {input.R2} 2> {log.err}
        """