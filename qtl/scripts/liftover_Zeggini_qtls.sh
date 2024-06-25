#!/bin/bash
#SBATCH -t 12:00:00
#SBATCH --mem=4G

module load ucsctools

liftOver /work/users/n/e/nekramer/CQTL/qtl/data/Steinberg_2021/raw/eQTL_LowGradeCartilage_hg19.bed \
    /proj/phanstiel_lab/References/genomes/hg19/hg19ToHg38.over.chain.gz \
    /work/users/n/e/nekramer/CQTL/qtl/data/Steinberg_2021/raw/eQTL_LowGradeCartilage_hg38.bed \
    /work/users/n/e/nekramer/CQTL/qtl/data/Steinberg_2021/raw/eQTL_LowGradeCartilage_hg19ToHg38.unmapped.bed \
    -bedPlus=1


liftOver /work/users/n/e/nekramer/CQTL/qtl/data/Steinberg_2021/raw/eQTL_HighGradeCartilage_hg19.bed \
    /proj/phanstiel_lab/References/genomes/hg19/hg19ToHg38.over.chain.gz \
    /work/users/n/e/nekramer/CQTL/qtl/data/Steinberg_2021/raw/eQTL_HighGradeCartilage_hg38.bed \
    /work/users/n/e/nekramer/CQTL/qtl/data/Steinberg_2021/raw/eQTL_HighGradeCartilage_hg19ToHg38.unmapped.bed \
    -bedPlus=1