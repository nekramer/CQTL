#!/usr/bin/env python3


with open('/work/users/n/e/nekramer/Data/CQTL/AI/freeze/AIsigCTL05variantIDs.txt') as file:
    CTL_snps = file.read().splitlines()

with open('/work/users/n/e/nekramer/Data/CQTL/AI/freeze/AIsigFNF05variantIDs.txt') as file:
    FNF_snps = file.read().splitlines()    


rule all:
    input:
        expand('/work/users/n/e/nekramer/Data/CQTL/AI/freeze/LD_buddies/{snp}_{condition}_0.8.ld', snp = CTL_snps, condition = "CTL"),
        expand('/work/users/n/e/nekramer/Data/CQTL/AI/freeze/LD_buddies/{snp}_{condition}_0.8.ld', snp = FNF_snps, condition = "FNF"),
        expand('data/{snp}_{condition}_ld_rsids.csv', snp = CTL_snps, condition = "CTL"),
        expand('data/{snp}_{condition}_ld_rsids.csv', snp = FNF_snps, condition = "FNF"),
        expand('data/{snp}_{condition}_ldbuddies_final.csv', snp = CTL_snps, condition = "CTL"),
        expand('data/{snp}_{condition}_ldbuddies_final.csv', snp = FNF_snps, condition = "FNF"),
        expand('data/{condition}_LDbuddies.csv', condition = ['CTL', 'FNF'])


rule getLD_buddies:
    input:
        snps = "/work/users/n/e/nekramer/Data/CQTL/AI/freeze/AIsig{condition}05variantIDs.txt"
    output:
        '/work/users/n/e/nekramer/Data/CQTL/AI/freeze/LD_buddies/{snp}_{condition}_0.8.ld'
    log:
        out = "logs/getLD_buddies_{condition}_{snp}.out",
        err = "logs/getLD_buddies_{condition}_{snp}.err"
    params:
        LDref = "/work/users/n/e/nekramer/Data/CQTL/geno/CQTL_COA_01_GDA8_COA2_01_COA3_01_GDA8_COA4_COA5_COA6_COA7_ALL_qc"
    shell:
        """
        module load plink
        for snp in `cat {input.snps}`
        do
            plink --bfile {params.LDref} --ld-snp ${{snp}} --ld-window 200000 --ld-window-kb 1000 --ld-window-r2 0.8 --r2 --out /work/users/n/e/nekramer/Data/CQTL/AI/freeze/LD_buddies/${{snp}}_{wildcards.condition}_0.8
        done
        """

rule getLDbuddy_rsIDs:
    input:
        rules.getLD_buddies.output
    output:
        'data/{snp}_{condition}_ld_rsids.csv'
    log:
        out = "logs/getLDbuddy_rsIDs_{condition}_{snp}.out",
        err = "logs/getLDbuddy_rsIDs_{condition}_{snp}.err"
    shell:
        """
        module load python/3.9.6
        python3 scripts/AI_LD_buddies_RSIDs.py /work/users/n/e/nekramer/Data/CQTL/AI/freeze/LD_buddies/{wildcards.snp}_{wildcards.condition}_0.8.ld {wildcards.snp} {wildcards.condition}
        """

rule join_LDbuddies:
    input:
        positions = rules.getLD_buddies.output,
        rsids = rules.getLDbuddy_rsIDs.output
    output:
        'data/{snp}_{condition}_ldbuddies_final.csv'
    log:
        out = "logs/join_LDbuddies_{condition}_{snp}.out",
        err = "logs/join_LDbuddies_{condition}_{snp}.err"
    shell:
        """
        module load r/4.1.3
        Rscript scripts/join_LDbuddies.R {input.rsids} {input.positions} {wildcards.snp} {wildcards.condition}
        """

rule join_conditionLDbuddies:
    input:
        lambda wildcards: expand('data/{snp}_{condition}_ldbuddies_final.csv', snp = eval(str(wildcards.condition) + '_snps'), condition = wildcards.condition)
    output:
        'data/{condition}_LDbuddies.csv'
    log:
        out = "logs/join_conditionLDbuddies_{condition}.out",
        err = "logs/join_conditionLDbuddies_{condition}.err"
    shell:
        """
        module load r/4.1.3
        Rscript scripts/join_allLD_buddies.R {wildcards.condition} {input}
        """