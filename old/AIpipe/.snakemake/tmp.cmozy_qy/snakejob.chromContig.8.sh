#!/bin/sh
# properties = {"type": "single", "rule": "chromContig", "local": false, "input": ["/proj/phanstiel_lab/Data/processed/CQTL/geno/CQTL_COA_01_GDA8_COA2_01_COA3_01_GDA8_COA4_COA5/vcf/CQTL_COA_01_GDA8_COA2_01_COA3_01_GDA8_COA4_COA5_ALL_qc.vcf.gz"], "output": ["output/vcf/CQTL_COA_01_GDA8_COA2_01_COA3_01_GDA8_COA4_COA5_renameContig.vcf"], "wildcards": {}, "params": {"chromNames": "chromNames.txt"}, "log": ["output/vcf/logs/renameContig.err"], "threads": 1, "resources": {}, "jobid": 8, "cluster": {"name": "chromContig", "partition": "general", "time": 4320, "cpusPerTask": "1", "memPerCpu": "4G", "nodes": 1, "output": "output/logs_slurm/chromContig.8.out", "error": "output/logs_slurm/chromContig.8.err"}}
 cd /pine/scr/n/e/nekramer/8_1_22_AIproc/CQTL/AIpipe && \
/pine/scr/n/e/nekramer/8_1_22_AIproc/CQTL/AIpipe/env/bin/python3 \
-m snakemake output/vcf/CQTL_COA_01_GDA8_COA2_01_COA3_01_GDA8_COA4_COA5_renameContig.vcf --snakefile /pine/scr/n/e/nekramer/8_1_22_AIproc/CQTL/AIpipe/snakefiles/VCFproc \
--force -j --keep-target-files --keep-remote \
--wait-for-files /pine/scr/n/e/nekramer/8_1_22_AIproc/CQTL/AIpipe/.snakemake/tmp.cmozy_qy /proj/phanstiel_lab/Data/processed/CQTL/geno/CQTL_COA_01_GDA8_COA2_01_COA3_01_GDA8_COA4_COA5/vcf/CQTL_COA_01_GDA8_COA2_01_COA3_01_GDA8_COA4_COA5_ALL_qc.vcf.gz --latency-wait 5 \
 --attempt 1 --force-use-threads \
--wrapper-prefix None \
--directory /pine/scr/n/e/nekramer/8_1_22_AIproc/CQTL/AIpipe --configfiles /pine/scr/n/e/nekramer/8_1_22_AIproc/CQTL/AIpipe/config/config.yaml  --allowed-rules chromContig --nocolor --notemp --no-hooks --nolock \
--mode 2  && exit 0 || exit 1

