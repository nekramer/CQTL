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
        zip = expand('output/{{group}}/qc/{{group}}_{R}_fastqc.zip', R=['R1', 'R2']),
        html = expand('output/{{group}}/qc/{{group}}_{R}_fastqc.html',R=['R1', 'R2'])
    threads: 2
    log:
        err = "output/{group}/logs/{group}_qc.err"
    shell:
        """
        module load fastqc/0.11.8
        mkdir -p output/{wildcards.group}/qc
        fastqc -t {threads} -o output/{wildcards.group}/qc {input.R1} {input.R2} 2> {log.err}
        """

rule trim:
    input:
        R1 = lambda wildcards: ['output/{group}/fastq/{group}_R1.fastq.gz'.format(group=wildcards.group)],
        R2 = lambda wildcards: ['output/{group}/fastq/{group}_R2.fastq.gz'.format(group=wildcards.group)]
    output:
        trim1 = temp("output/{group}/trim/{group}_R1_val_1.fq.gz"),
        trim2 = temp("output/{group}/trim/{group}_R2_val_2.fq.gz")
    threads: 4
    log:
        err = "output/{group}/logs/{group}_trim.err"
    shell:
        """
        module load trim_galore/0.6.2
        module load python/3.6.6
        module load pigz
        mkdir -p output/{wildcards.group}/trim
        trim_galore -o output/{wildcards.group}/trim --cores {threads} --path_to_cutadapt /nas/longleaf/apps/cutadapt/2.9/venv/bin/cutadapt --paired {input.R1} {input.R2} 2> {log.err}
        """