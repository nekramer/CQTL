#!/bin/bash

module load r/4.2.1

for chr in {1..22}
do
    awk '{print $1 " " $2 " " $2 " " $3 " " $4}' genetic_map_GRCh37_chr${chr}.txt > genetic_map_GRCh37_chr${chr}.bed
    sed -i '1d' genetic_map_GRCh37_chr${chr}.bed
    /proj/phanstiel_lab/Software/liftOver genetic_map_GRCh37_chr${chr}.bed hg19ToHg38.over.chain genetic_map_GRCh38_chr${chr}.bed genetic_map_GRCh38_chr${chr}.unMapped
    awk '{print $2 " " $4 " " $5}' genetic_map_GRCh38_chr${chr}.bed > genetic_map_GRCh38_chr${chr}.txt
    sed -i '1s/^/PhysicalPosition(bp) RecombinationRate GenomicPosition(cM)\n/' genetic_map_GRCh38_chr${chr}.txt

    Rscript sortLiftedFile.R ${chr}
    
    
done