#!/bin/sh
# properties = {"type": "single", "rule": "mergeAlleleCounts", "local": false, "input": ["output/CQTL_AM7180_R_CTL_18_1/alleleCounts/CQTL_AM7180_R_CTL_18_1_alleleCounts.csv", "output/CQTL_AM7180_R_FNF_18_1/alleleCounts/CQTL_AM7180_R_FNF_18_1_alleleCounts.csv", "output/CQTL_AM7181_R_CTL_18_1/alleleCounts/CQTL_AM7181_R_CTL_18_1_alleleCounts.csv", "output/CQTL_AM7181_R_FNF_18_1/alleleCounts/CQTL_AM7181_R_FNF_18_1_alleleCounts.csv", "output/CQTL_AM7182_R_CTL_18_1/alleleCounts/CQTL_AM7182_R_CTL_18_1_alleleCounts.csv", "output/CQTL_AM7182_R_FNF_18_1/alleleCounts/CQTL_AM7182_R_FNF_18_1_alleleCounts.csv", "output/CQTL_AM7188_R_CTL_18_1/alleleCounts/CQTL_AM7188_R_CTL_18_1_alleleCounts.csv", "output/CQTL_AM7188_R_FNF_18_1/alleleCounts/CQTL_AM7188_R_FNF_18_1_alleleCounts.csv", "output/CQTL_AM7189_R_CTL_18_1/alleleCounts/CQTL_AM7189_R_CTL_18_1_alleleCounts.csv", "output/CQTL_AM7189_R_FNF_18_1/alleleCounts/CQTL_AM7189_R_FNF_18_1_alleleCounts.csv", "output/CQTL_AM7196_R_CTL_18_1/alleleCounts/CQTL_AM7196_R_CTL_18_1_alleleCounts.csv", "output/CQTL_AM7196_R_FNF_18_1/alleleCounts/CQTL_AM7196_R_FNF_18_1_alleleCounts.csv", "output/CQTL_AM7197_R_CTL_18_1/alleleCounts/CQTL_AM7197_R_CTL_18_1_alleleCounts.csv", "output/CQTL_AM7197_R_FNF_18_1/alleleCounts/CQTL_AM7197_R_FNF_18_1_alleleCounts.csv", "output/CQTL_AM7203_R_CTL_18_1/alleleCounts/CQTL_AM7203_R_CTL_18_1_alleleCounts.csv", "output/CQTL_AM7203_R_FNF_18_1/alleleCounts/CQTL_AM7203_R_FNF_18_1_alleleCounts.csv", "output/CQTL_AM7204_R_CTL_18_1/alleleCounts/CQTL_AM7204_R_CTL_18_1_alleleCounts.csv", "output/CQTL_AM7204_R_FNF_18_1/alleleCounts/CQTL_AM7204_R_FNF_18_1_alleleCounts.csv", "output/CQTL_AM7205_R_CTL_18_1/alleleCounts/CQTL_AM7205_R_CTL_18_1_alleleCounts.csv", "output/CQTL_AM7205_R_FNF_18_1/alleleCounts/CQTL_AM7205_R_FNF_18_1_alleleCounts.csv", "output/CQTL_AM7208_R_CTL_18_1/alleleCounts/CQTL_AM7208_R_CTL_18_1_alleleCounts.csv", "output/CQTL_AM7208_R_FNF_18_1/alleleCounts/CQTL_AM7208_R_FNF_18_1_alleleCounts.csv", "output/CQTL_AM7209_R_CTL_18_1/alleleCounts/CQTL_AM7209_R_CTL_18_1_alleleCounts.csv", "output/CQTL_AM7209_R_FNF_18_1/alleleCounts/CQTL_AM7209_R_FNF_18_1_alleleCounts.csv", "output/CQTL_AM7211_R_CTL_18_1/alleleCounts/CQTL_AM7211_R_CTL_18_1_alleleCounts.csv", "output/CQTL_AM7211_R_FNF_18_1/alleleCounts/CQTL_AM7211_R_FNF_18_1_alleleCounts.csv", "output/CQTL_AM7213_R_CTL_18_1/alleleCounts/CQTL_AM7213_R_CTL_18_1_alleleCounts.csv", "output/CQTL_AM7213_R_FNF_18_1/alleleCounts/CQTL_AM7213_R_FNF_18_1_alleleCounts.csv", "output/CQTL_AM7214_R_CTL_18_1/alleleCounts/CQTL_AM7214_R_CTL_18_1_alleleCounts.csv", "output/CQTL_AM7214_R_FNF_18_1/alleleCounts/CQTL_AM7214_R_FNF_18_1_alleleCounts.csv", "output/CQTL_AM7215_R_CTL_18_1/alleleCounts/CQTL_AM7215_R_CTL_18_1_alleleCounts.csv", "output/CQTL_AM7215_R_FNF_18_1/alleleCounts/CQTL_AM7215_R_FNF_18_1_alleleCounts.csv", "output/CQTL_AM7216_R_CTL_18_1/alleleCounts/CQTL_AM7216_R_CTL_18_1_alleleCounts.csv", "output/CQTL_AM7216_R_FNF_18_1/alleleCounts/CQTL_AM7216_R_FNF_18_1_alleleCounts.csv", "output/CQTL_AM7221_R_CTL_18_1/alleleCounts/CQTL_AM7221_R_CTL_18_1_alleleCounts.csv", "output/CQTL_AM7221_R_FNF_18_1/alleleCounts/CQTL_AM7221_R_FNF_18_1_alleleCounts.csv", "output/CQTL_AM7223_R_CTL_18_1/alleleCounts/CQTL_AM7223_R_CTL_18_1_alleleCounts.csv", "output/CQTL_AM7223_R_FNF_18_1/alleleCounts/CQTL_AM7223_R_FNF_18_1_alleleCounts.csv", "output/CQTL_AM7224_R_CTL_18_1/alleleCounts/CQTL_AM7224_R_CTL_18_1_alleleCounts.csv", "output/CQTL_AM7224_R_FNF_18_1/alleleCounts/CQTL_AM7224_R_FNF_18_1_alleleCounts.csv", "output/CQTL_AM7226_R_CTL_18_1/alleleCounts/CQTL_AM7226_R_CTL_18_1_alleleCounts.csv", "output/CQTL_AM7226_R_FNF_18_1/alleleCounts/CQTL_AM7226_R_FNF_18_1_alleleCounts.csv", "output/CQTL_AM7228_R_CTL_18_1/alleleCounts/CQTL_AM7228_R_CTL_18_1_alleleCounts.csv", "output/CQTL_AM7228_R_FNF_18_1/alleleCounts/CQTL_AM7228_R_FNF_18_1_alleleCounts.csv", "output/CQTL_AM7229_R_CTL_18_1/alleleCounts/CQTL_AM7229_R_CTL_18_1_alleleCounts.csv", "output/CQTL_AM7229_R_FNF_18_1/alleleCounts/CQTL_AM7229_R_FNF_18_1_alleleCounts.csv", "output/CQTL_AM7230_R_CTL_18_1/alleleCounts/CQTL_AM7230_R_CTL_18_1_alleleCounts.csv", "output/CQTL_AM7230_R_FNF_18_1/alleleCounts/CQTL_AM7230_R_FNF_18_1_alleleCounts.csv", "output/CQTL_AM7236_R_CTL_18_1/alleleCounts/CQTL_AM7236_R_CTL_18_1_alleleCounts.csv", "output/CQTL_AM7236_R_FNF_18_1/alleleCounts/CQTL_AM7236_R_FNF_18_1_alleleCounts.csv", "output/CQTL_AM7237_R_CTL_18_1/alleleCounts/CQTL_AM7237_R_CTL_18_1_alleleCounts.csv", "output/CQTL_AM7237_R_FNF_18_1/alleleCounts/CQTL_AM7237_R_FNF_18_1_alleleCounts.csv", "output/CQTL_AM7241_R_CTL_18_1/alleleCounts/CQTL_AM7241_R_CTL_18_1_alleleCounts.csv", "output/CQTL_AM7241_R_FNF_18_1/alleleCounts/CQTL_AM7241_R_FNF_18_1_alleleCounts.csv", "output/CQTL_AM7242_R_CTL_18_1/alleleCounts/CQTL_AM7242_R_CTL_18_1_alleleCounts.csv", "output/CQTL_AM7242_R_FNF_18_1/alleleCounts/CQTL_AM7242_R_FNF_18_1_alleleCounts.csv", "output/CQTL_AM7243_R_CTL_18_1/alleleCounts/CQTL_AM7243_R_CTL_18_1_alleleCounts.csv", "output/CQTL_AM7243_R_FNF_18_1/alleleCounts/CQTL_AM7243_R_FNF_18_1_alleleCounts.csv", "output/CQTL_AM7244_R_CTL_18_1/alleleCounts/CQTL_AM7244_R_CTL_18_1_alleleCounts.csv", "output/CQTL_AM7244_R_FNF_18_1/alleleCounts/CQTL_AM7244_R_FNF_18_1_alleleCounts.csv", "output/CQTL_AM7255_R_CTL_18_1/alleleCounts/CQTL_AM7255_R_CTL_18_1_alleleCounts.csv", "output/CQTL_AM7255_R_FNF_18_1/alleleCounts/CQTL_AM7255_R_FNF_18_1_alleleCounts.csv", "output/CQTL_AM7256_R_CTL_18_1/alleleCounts/CQTL_AM7256_R_CTL_18_1_alleleCounts.csv", "output/CQTL_AM7256_R_FNF_18_1/alleleCounts/CQTL_AM7256_R_FNF_18_1_alleleCounts.csv", "output/CQTL_AM7260_R_CTL_18_1/alleleCounts/CQTL_AM7260_R_CTL_18_1_alleleCounts.csv", "output/CQTL_AM7260_R_FNF_18_1/alleleCounts/CQTL_AM7260_R_FNF_18_1_alleleCounts.csv", "output/CQTL_AM7261_R_CTL_18_1/alleleCounts/CQTL_AM7261_R_CTL_18_1_alleleCounts.csv", "output/CQTL_AM7261_R_FNF_18_1/alleleCounts/CQTL_AM7261_R_FNF_18_1_alleleCounts.csv", "output/CQTL_AM7266_R_CTL_18_1/alleleCounts/CQTL_AM7266_R_CTL_18_1_alleleCounts.csv", "output/CQTL_AM7266_R_FNF_18_1/alleleCounts/CQTL_AM7266_R_FNF_18_1_alleleCounts.csv", "output/CQTL_AM7269_R_CTL_18_1/alleleCounts/CQTL_AM7269_R_CTL_18_1_alleleCounts.csv", "output/CQTL_AM7269_R_FNF_18_1/alleleCounts/CQTL_AM7269_R_FNF_18_1_alleleCounts.csv", "output/CQTL_AM7270_R_CTL_18_1/alleleCounts/CQTL_AM7270_R_CTL_18_1_alleleCounts.csv", "output/CQTL_AM7270_R_FNF_18_1/alleleCounts/CQTL_AM7270_R_FNF_18_1_alleleCounts.csv", "output/CQTL_AM7272_R_CTL_18_1/alleleCounts/CQTL_AM7272_R_CTL_18_1_alleleCounts.csv", "output/CQTL_AM7272_R_FNF_18_1/alleleCounts/CQTL_AM7272_R_FNF_18_1_alleleCounts.csv", "output/CQTL_AM7273_R_CTL_18_2/alleleCounts/CQTL_AM7273_R_CTL_18_2_alleleCounts.csv", "output/CQTL_AM7273_R_FNF_18_1/alleleCounts/CQTL_AM7273_R_FNF_18_1_alleleCounts.csv", "output/CQTL_AM7277_R_CTL_18_1/alleleCounts/CQTL_AM7277_R_CTL_18_1_alleleCounts.csv", "output/CQTL_AM7277_R_FNF_18_1/alleleCounts/CQTL_AM7277_R_FNF_18_1_alleleCounts.csv", "output/CQTL_AM7278_R_CTL_18_1/alleleCounts/CQTL_AM7278_R_CTL_18_1_alleleCounts.csv", "output/CQTL_AM7278_R_FNF_18_1/alleleCounts/CQTL_AM7278_R_FNF_18_1_alleleCounts.csv", "output/CQTL_AM7280_R_CTL_18_1/alleleCounts/CQTL_AM7280_R_CTL_18_1_alleleCounts.csv", "output/CQTL_AM7280_R_FNF_18_1/alleleCounts/CQTL_AM7280_R_FNF_18_1_alleleCounts.csv", "output/CQTL_AM7283_R_CTL_18_1/alleleCounts/CQTL_AM7283_R_CTL_18_1_alleleCounts.csv", "output/CQTL_AM7283_R_FNF_18_1/alleleCounts/CQTL_AM7283_R_FNF_18_1_alleleCounts.csv", "output/CQTL_AM7284_R_CTL_18_1/alleleCounts/CQTL_AM7284_R_CTL_18_1_alleleCounts.csv", "output/CQTL_AM7284_R_FNF_18_1/alleleCounts/CQTL_AM7284_R_FNF_18_1_alleleCounts.csv", "output/CQTL_AM7285_R_CTL_18_1/alleleCounts/CQTL_AM7285_R_CTL_18_1_alleleCounts.csv", "output/CQTL_AM7285_R_FNF_18_1/alleleCounts/CQTL_AM7285_R_FNF_18_1_alleleCounts.csv", "output/CQTL_AM7294_R_CTL_18_1/alleleCounts/CQTL_AM7294_R_CTL_18_1_alleleCounts.csv", "output/CQTL_AM7294_R_FNF_18_1/alleleCounts/CQTL_AM7294_R_FNF_18_1_alleleCounts.csv", "output/CQTL_AM7295_R_CTL_18_1/alleleCounts/CQTL_AM7295_R_CTL_18_1_alleleCounts.csv", "output/CQTL_AM7295_R_FNF_18_1/alleleCounts/CQTL_AM7295_R_FNF_18_1_alleleCounts.csv", "output/CQTL_AM7302_R_CTL_18_1/alleleCounts/CQTL_AM7302_R_CTL_18_1_alleleCounts.csv", "output/CQTL_AM7302_R_FNF_18_1/alleleCounts/CQTL_AM7302_R_FNF_18_1_alleleCounts.csv", "output/CQTL_AM7303_R_CTL_18_1/alleleCounts/CQTL_AM7303_R_CTL_18_1_alleleCounts.csv", "output/CQTL_AM7303_R_FNF_18_1/alleleCounts/CQTL_AM7303_R_FNF_18_1_alleleCounts.csv", "output/CQTL_AM7304_R_CTL_18_1/alleleCounts/CQTL_AM7304_R_CTL_18_1_alleleCounts.csv", "output/CQTL_AM7304_R_FNF_18_1/alleleCounts/CQTL_AM7304_R_FNF_18_1_alleleCounts.csv", "output/CQTL_AM7312_R_CTL_18_1/alleleCounts/CQTL_AM7312_R_CTL_18_1_alleleCounts.csv", "output/CQTL_AM7312_R_FNF_18_1/alleleCounts/CQTL_AM7312_R_FNF_18_1_alleleCounts.csv", "output/CQTL_AM7313_R_CTL_18_1/alleleCounts/CQTL_AM7313_R_CTL_18_1_alleleCounts.csv", "output/CQTL_AM7313_R_FNF_18_1/alleleCounts/CQTL_AM7313_R_FNF_18_1_alleleCounts.csv", "output/CQTL_AM7318_R_CTL_18_1/alleleCounts/CQTL_AM7318_R_CTL_18_1_alleleCounts.csv", "output/CQTL_AM7318_R_FNF_18_1/alleleCounts/CQTL_AM7318_R_FNF_18_1_alleleCounts.csv", "output/CQTL_AM7319_R_CTL_18_1/alleleCounts/CQTL_AM7319_R_CTL_18_1_alleleCounts.csv", "output/CQTL_AM7319_R_FNF_18_1/alleleCounts/CQTL_AM7319_R_FNF_18_1_alleleCounts.csv", "output/CQTL_AM7320_R_CTL_18_1/alleleCounts/CQTL_AM7320_R_CTL_18_1_alleleCounts.csv", "output/CQTL_AM7320_R_FNF_18_1/alleleCounts/CQTL_AM7320_R_FNF_18_1_alleleCounts.csv", "output/CQTL_AM7321_R_CTL_18_1/alleleCounts/CQTL_AM7321_R_CTL_18_1_alleleCounts.csv", "output/CQTL_AM7321_R_FNF_18_1/alleleCounts/CQTL_AM7321_R_FNF_18_1_alleleCounts.csv"], "output": ["output/AI/alleleCountsMatrix.txt", "output/AI/weightsMatrix.txt", "output/AI/colData.txt"], "wildcards": {}, "params": {"donors": "AM7180,AM7181,AM7182,AM7188,AM7189,AM7196,AM7197,AM7203,AM7204,AM7205,AM7208,AM7209,AM7211,AM7213,AM7214,AM7215,AM7216,AM7221,AM7223,AM7224,AM7226,AM7228,AM7229,AM7230,AM7236,AM7237,AM7241,AM7242,AM7243,AM7244,AM7255,AM7256,AM7260,AM7261,AM7266,AM7269,AM7270,AM7272,AM7273,AM7277,AM7278,AM7280,AM7283,AM7284,AM7285,AM7294,AM7295,AM7302,AM7303,AM7304,AM7312,AM7313,AM7318,AM7319,AM7320,AM7321", "conditions": "CTL,FNF"}, "log": ["output/AI/logs/concatAlleleCounts.out"], "threads": 1, "resources": {}, "jobid": 1, "cluster": {"name": "mergeAlleleCounts", "partition": "general", "time": 4320, "cpusPerTask": "1", "memPerCpu": "8G", "nodes": 1, "output": "output/logs_slurm/mergeAlleleCounts.1.out", "error": "output/logs_slurm/mergeAlleleCounts.1.err"}}
 cd /pine/scr/n/e/nekramer/8_1_22_AIproc/CQTL/AIpipe && \
/pine/scr/n/e/nekramer/8_1_22_AIproc/CQTL/AIpipe/env/bin/python3 \
-m snakemake output/AI/alleleCountsMatrix.txt --snakefile /pine/scr/n/e/nekramer/8_1_22_AIproc/CQTL/AIpipe/snakefiles/NEW_AIanalysis.smk \
--force -j --keep-target-files --keep-remote \
--wait-for-files /pine/scr/n/e/nekramer/8_1_22_AIproc/CQTL/AIpipe/.snakemake/tmp.muqjnd8o output/CQTL_AM7180_R_CTL_18_1/alleleCounts/CQTL_AM7180_R_CTL_18_1_alleleCounts.csv output/CQTL_AM7180_R_FNF_18_1/alleleCounts/CQTL_AM7180_R_FNF_18_1_alleleCounts.csv output/CQTL_AM7181_R_CTL_18_1/alleleCounts/CQTL_AM7181_R_CTL_18_1_alleleCounts.csv output/CQTL_AM7181_R_FNF_18_1/alleleCounts/CQTL_AM7181_R_FNF_18_1_alleleCounts.csv output/CQTL_AM7182_R_CTL_18_1/alleleCounts/CQTL_AM7182_R_CTL_18_1_alleleCounts.csv output/CQTL_AM7182_R_FNF_18_1/alleleCounts/CQTL_AM7182_R_FNF_18_1_alleleCounts.csv output/CQTL_AM7188_R_CTL_18_1/alleleCounts/CQTL_AM7188_R_CTL_18_1_alleleCounts.csv output/CQTL_AM7188_R_FNF_18_1/alleleCounts/CQTL_AM7188_R_FNF_18_1_alleleCounts.csv output/CQTL_AM7189_R_CTL_18_1/alleleCounts/CQTL_AM7189_R_CTL_18_1_alleleCounts.csv output/CQTL_AM7189_R_FNF_18_1/alleleCounts/CQTL_AM7189_R_FNF_18_1_alleleCounts.csv output/CQTL_AM7196_R_CTL_18_1/alleleCounts/CQTL_AM7196_R_CTL_18_1_alleleCounts.csv output/CQTL_AM7196_R_FNF_18_1/alleleCounts/CQTL_AM7196_R_FNF_18_1_alleleCounts.csv output/CQTL_AM7197_R_CTL_18_1/alleleCounts/CQTL_AM7197_R_CTL_18_1_alleleCounts.csv output/CQTL_AM7197_R_FNF_18_1/alleleCounts/CQTL_AM7197_R_FNF_18_1_alleleCounts.csv output/CQTL_AM7203_R_CTL_18_1/alleleCounts/CQTL_AM7203_R_CTL_18_1_alleleCounts.csv output/CQTL_AM7203_R_FNF_18_1/alleleCounts/CQTL_AM7203_R_FNF_18_1_alleleCounts.csv output/CQTL_AM7204_R_CTL_18_1/alleleCounts/CQTL_AM7204_R_CTL_18_1_alleleCounts.csv output/CQTL_AM7204_R_FNF_18_1/alleleCounts/CQTL_AM7204_R_FNF_18_1_alleleCounts.csv output/CQTL_AM7205_R_CTL_18_1/alleleCounts/CQTL_AM7205_R_CTL_18_1_alleleCounts.csv output/CQTL_AM7205_R_FNF_18_1/alleleCounts/CQTL_AM7205_R_FNF_18_1_alleleCounts.csv output/CQTL_AM7208_R_CTL_18_1/alleleCounts/CQTL_AM7208_R_CTL_18_1_alleleCounts.csv output/CQTL_AM7208_R_FNF_18_1/alleleCounts/CQTL_AM7208_R_FNF_18_1_alleleCounts.csv output/CQTL_AM7209_R_CTL_18_1/alleleCounts/CQTL_AM7209_R_CTL_18_1_alleleCounts.csv output/CQTL_AM7209_R_FNF_18_1/alleleCounts/CQTL_AM7209_R_FNF_18_1_alleleCounts.csv output/CQTL_AM7211_R_CTL_18_1/alleleCounts/CQTL_AM7211_R_CTL_18_1_alleleCounts.csv output/CQTL_AM7211_R_FNF_18_1/alleleCounts/CQTL_AM7211_R_FNF_18_1_alleleCounts.csv output/CQTL_AM7213_R_CTL_18_1/alleleCounts/CQTL_AM7213_R_CTL_18_1_alleleCounts.csv output/CQTL_AM7213_R_FNF_18_1/alleleCounts/CQTL_AM7213_R_FNF_18_1_alleleCounts.csv output/CQTL_AM7214_R_CTL_18_1/alleleCounts/CQTL_AM7214_R_CTL_18_1_alleleCounts.csv output/CQTL_AM7214_R_FNF_18_1/alleleCounts/CQTL_AM7214_R_FNF_18_1_alleleCounts.csv output/CQTL_AM7215_R_CTL_18_1/alleleCounts/CQTL_AM7215_R_CTL_18_1_alleleCounts.csv output/CQTL_AM7215_R_FNF_18_1/alleleCounts/CQTL_AM7215_R_FNF_18_1_alleleCounts.csv output/CQTL_AM7216_R_CTL_18_1/alleleCounts/CQTL_AM7216_R_CTL_18_1_alleleCounts.csv output/CQTL_AM7216_R_FNF_18_1/alleleCounts/CQTL_AM7216_R_FNF_18_1_alleleCounts.csv output/CQTL_AM7221_R_CTL_18_1/alleleCounts/CQTL_AM7221_R_CTL_18_1_alleleCounts.csv output/CQTL_AM7221_R_FNF_18_1/alleleCounts/CQTL_AM7221_R_FNF_18_1_alleleCounts.csv output/CQTL_AM7223_R_CTL_18_1/alleleCounts/CQTL_AM7223_R_CTL_18_1_alleleCounts.csv output/CQTL_AM7223_R_FNF_18_1/alleleCounts/CQTL_AM7223_R_FNF_18_1_alleleCounts.csv output/CQTL_AM7224_R_CTL_18_1/alleleCounts/CQTL_AM7224_R_CTL_18_1_alleleCounts.csv output/CQTL_AM7224_R_FNF_18_1/alleleCounts/CQTL_AM7224_R_FNF_18_1_alleleCounts.csv output/CQTL_AM7226_R_CTL_18_1/alleleCounts/CQTL_AM7226_R_CTL_18_1_alleleCounts.csv output/CQTL_AM7226_R_FNF_18_1/alleleCounts/CQTL_AM7226_R_FNF_18_1_alleleCounts.csv output/CQTL_AM7228_R_CTL_18_1/alleleCounts/CQTL_AM7228_R_CTL_18_1_alleleCounts.csv output/CQTL_AM7228_R_FNF_18_1/alleleCounts/CQTL_AM7228_R_FNF_18_1_alleleCounts.csv output/CQTL_AM7229_R_CTL_18_1/alleleCounts/CQTL_AM7229_R_CTL_18_1_alleleCounts.csv output/CQTL_AM7229_R_FNF_18_1/alleleCounts/CQTL_AM7229_R_FNF_18_1_alleleCounts.csv output/CQTL_AM7230_R_CTL_18_1/alleleCounts/CQTL_AM7230_R_CTL_18_1_alleleCounts.csv output/CQTL_AM7230_R_FNF_18_1/alleleCounts/CQTL_AM7230_R_FNF_18_1_alleleCounts.csv output/CQTL_AM7236_R_CTL_18_1/alleleCounts/CQTL_AM7236_R_CTL_18_1_alleleCounts.csv output/CQTL_AM7236_R_FNF_18_1/alleleCounts/CQTL_AM7236_R_FNF_18_1_alleleCounts.csv output/CQTL_AM7237_R_CTL_18_1/alleleCounts/CQTL_AM7237_R_CTL_18_1_alleleCounts.csv output/CQTL_AM7237_R_FNF_18_1/alleleCounts/CQTL_AM7237_R_FNF_18_1_alleleCounts.csv output/CQTL_AM7241_R_CTL_18_1/alleleCounts/CQTL_AM7241_R_CTL_18_1_alleleCounts.csv output/CQTL_AM7241_R_FNF_18_1/alleleCounts/CQTL_AM7241_R_FNF_18_1_alleleCounts.csv output/CQTL_AM7242_R_CTL_18_1/alleleCounts/CQTL_AM7242_R_CTL_18_1_alleleCounts.csv output/CQTL_AM7242_R_FNF_18_1/alleleCounts/CQTL_AM7242_R_FNF_18_1_alleleCounts.csv output/CQTL_AM7243_R_CTL_18_1/alleleCounts/CQTL_AM7243_R_CTL_18_1_alleleCounts.csv output/CQTL_AM7243_R_FNF_18_1/alleleCounts/CQTL_AM7243_R_FNF_18_1_alleleCounts.csv output/CQTL_AM7244_R_CTL_18_1/alleleCounts/CQTL_AM7244_R_CTL_18_1_alleleCounts.csv output/CQTL_AM7244_R_FNF_18_1/alleleCounts/CQTL_AM7244_R_FNF_18_1_alleleCounts.csv output/CQTL_AM7255_R_CTL_18_1/alleleCounts/CQTL_AM7255_R_CTL_18_1_alleleCounts.csv output/CQTL_AM7255_R_FNF_18_1/alleleCounts/CQTL_AM7255_R_FNF_18_1_alleleCounts.csv output/CQTL_AM7256_R_CTL_18_1/alleleCounts/CQTL_AM7256_R_CTL_18_1_alleleCounts.csv output/CQTL_AM7256_R_FNF_18_1/alleleCounts/CQTL_AM7256_R_FNF_18_1_alleleCounts.csv output/CQTL_AM7260_R_CTL_18_1/alleleCounts/CQTL_AM7260_R_CTL_18_1_alleleCounts.csv output/CQTL_AM7260_R_FNF_18_1/alleleCounts/CQTL_AM7260_R_FNF_18_1_alleleCounts.csv output/CQTL_AM7261_R_CTL_18_1/alleleCounts/CQTL_AM7261_R_CTL_18_1_alleleCounts.csv output/CQTL_AM7261_R_FNF_18_1/alleleCounts/CQTL_AM7261_R_FNF_18_1_alleleCounts.csv output/CQTL_AM7266_R_CTL_18_1/alleleCounts/CQTL_AM7266_R_CTL_18_1_alleleCounts.csv output/CQTL_AM7266_R_FNF_18_1/alleleCounts/CQTL_AM7266_R_FNF_18_1_alleleCounts.csv output/CQTL_AM7269_R_CTL_18_1/alleleCounts/CQTL_AM7269_R_CTL_18_1_alleleCounts.csv output/CQTL_AM7269_R_FNF_18_1/alleleCounts/CQTL_AM7269_R_FNF_18_1_alleleCounts.csv output/CQTL_AM7270_R_CTL_18_1/alleleCounts/CQTL_AM7270_R_CTL_18_1_alleleCounts.csv output/CQTL_AM7270_R_FNF_18_1/alleleCounts/CQTL_AM7270_R_FNF_18_1_alleleCounts.csv output/CQTL_AM7272_R_CTL_18_1/alleleCounts/CQTL_AM7272_R_CTL_18_1_alleleCounts.csv output/CQTL_AM7272_R_FNF_18_1/alleleCounts/CQTL_AM7272_R_FNF_18_1_alleleCounts.csv output/CQTL_AM7273_R_CTL_18_2/alleleCounts/CQTL_AM7273_R_CTL_18_2_alleleCounts.csv output/CQTL_AM7273_R_FNF_18_1/alleleCounts/CQTL_AM7273_R_FNF_18_1_alleleCounts.csv output/CQTL_AM7277_R_CTL_18_1/alleleCounts/CQTL_AM7277_R_CTL_18_1_alleleCounts.csv output/CQTL_AM7277_R_FNF_18_1/alleleCounts/CQTL_AM7277_R_FNF_18_1_alleleCounts.csv output/CQTL_AM7278_R_CTL_18_1/alleleCounts/CQTL_AM7278_R_CTL_18_1_alleleCounts.csv output/CQTL_AM7278_R_FNF_18_1/alleleCounts/CQTL_AM7278_R_FNF_18_1_alleleCounts.csv output/CQTL_AM7280_R_CTL_18_1/alleleCounts/CQTL_AM7280_R_CTL_18_1_alleleCounts.csv output/CQTL_AM7280_R_FNF_18_1/alleleCounts/CQTL_AM7280_R_FNF_18_1_alleleCounts.csv output/CQTL_AM7283_R_CTL_18_1/alleleCounts/CQTL_AM7283_R_CTL_18_1_alleleCounts.csv output/CQTL_AM7283_R_FNF_18_1/alleleCounts/CQTL_AM7283_R_FNF_18_1_alleleCounts.csv output/CQTL_AM7284_R_CTL_18_1/alleleCounts/CQTL_AM7284_R_CTL_18_1_alleleCounts.csv output/CQTL_AM7284_R_FNF_18_1/alleleCounts/CQTL_AM7284_R_FNF_18_1_alleleCounts.csv output/CQTL_AM7285_R_CTL_18_1/alleleCounts/CQTL_AM7285_R_CTL_18_1_alleleCounts.csv output/CQTL_AM7285_R_FNF_18_1/alleleCounts/CQTL_AM7285_R_FNF_18_1_alleleCounts.csv output/CQTL_AM7294_R_CTL_18_1/alleleCounts/CQTL_AM7294_R_CTL_18_1_alleleCounts.csv output/CQTL_AM7294_R_FNF_18_1/alleleCounts/CQTL_AM7294_R_FNF_18_1_alleleCounts.csv output/CQTL_AM7295_R_CTL_18_1/alleleCounts/CQTL_AM7295_R_CTL_18_1_alleleCounts.csv output/CQTL_AM7295_R_FNF_18_1/alleleCounts/CQTL_AM7295_R_FNF_18_1_alleleCounts.csv output/CQTL_AM7302_R_CTL_18_1/alleleCounts/CQTL_AM7302_R_CTL_18_1_alleleCounts.csv output/CQTL_AM7302_R_FNF_18_1/alleleCounts/CQTL_AM7302_R_FNF_18_1_alleleCounts.csv output/CQTL_AM7303_R_CTL_18_1/alleleCounts/CQTL_AM7303_R_CTL_18_1_alleleCounts.csv output/CQTL_AM7303_R_FNF_18_1/alleleCounts/CQTL_AM7303_R_FNF_18_1_alleleCounts.csv output/CQTL_AM7304_R_CTL_18_1/alleleCounts/CQTL_AM7304_R_CTL_18_1_alleleCounts.csv output/CQTL_AM7304_R_FNF_18_1/alleleCounts/CQTL_AM7304_R_FNF_18_1_alleleCounts.csv output/CQTL_AM7312_R_CTL_18_1/alleleCounts/CQTL_AM7312_R_CTL_18_1_alleleCounts.csv output/CQTL_AM7312_R_FNF_18_1/alleleCounts/CQTL_AM7312_R_FNF_18_1_alleleCounts.csv output/CQTL_AM7313_R_CTL_18_1/alleleCounts/CQTL_AM7313_R_CTL_18_1_alleleCounts.csv output/CQTL_AM7313_R_FNF_18_1/alleleCounts/CQTL_AM7313_R_FNF_18_1_alleleCounts.csv output/CQTL_AM7318_R_CTL_18_1/alleleCounts/CQTL_AM7318_R_CTL_18_1_alleleCounts.csv output/CQTL_AM7318_R_FNF_18_1/alleleCounts/CQTL_AM7318_R_FNF_18_1_alleleCounts.csv output/CQTL_AM7319_R_CTL_18_1/alleleCounts/CQTL_AM7319_R_CTL_18_1_alleleCounts.csv output/CQTL_AM7319_R_FNF_18_1/alleleCounts/CQTL_AM7319_R_FNF_18_1_alleleCounts.csv output/CQTL_AM7320_R_CTL_18_1/alleleCounts/CQTL_AM7320_R_CTL_18_1_alleleCounts.csv output/CQTL_AM7320_R_FNF_18_1/alleleCounts/CQTL_AM7320_R_FNF_18_1_alleleCounts.csv output/CQTL_AM7321_R_CTL_18_1/alleleCounts/CQTL_AM7321_R_CTL_18_1_alleleCounts.csv output/CQTL_AM7321_R_FNF_18_1/alleleCounts/CQTL_AM7321_R_FNF_18_1_alleleCounts.csv --latency-wait 5 \
 --attempt 1 --force-use-threads \
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
   --allowed-rules mergeAlleleCounts --nocolor --notemp --no-hooks --nolock \
--mode 2  && exit 0 || exit 1

