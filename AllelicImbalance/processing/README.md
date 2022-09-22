# Allelic Imbalance processing

Snakemake workflows for preparing imputed genotyping data and RNA-seq data for 
[GATK ASERead Counter](https://gatk.broadinstitute.org/hc/en-us/articles/360037428291-ASEReadCounter) and preparing data for 
differential analysis.
Before running this pipeline, refer to `Genopipe` to QC and impute raw
genotyping data.

## Workflow

1. Edit the comma-separated `samplesheet.csv` with the names of `Read1` and `Read2` gzipped fastq files and the path to these files
under the `Sequencing_Directory` column. Additional sample metadata includes: `Proj`, `Donor`, `Condition`, `Time`, `Tech_Rep`, and `Seq_Rep`.

2. Edit `config/config_AIprocessing.yaml` for analysis-specific parameters:

    ```yaml
    ## Path to RNA-seq sample samplesheet
    samplesheet: 'samplesheet.csv'

    ## Path to imputed genotyping vcf
    vcf: '/path/to/imputed/genotyping.vcf.gz'

    ## Genome-specific reference parameters
    chromNames: 'chromNames.txt' # 2-column text file. 1st column: format of chromosome names in vcf file, 2nd column: required format of chromosome names for compatibility with sequencing contigs.
    genomeDir: '/proj/seq/data/STAR_genomes/GRCh38_p10_GENCODE' # Path to folder of STAR genome for alignment.
    chromSizes: '/proj/seq/data/STAR_genomes/GRCh38_p10_GENCODE/chrNameLength.txt' # Path to file of chromosome names and lengths.
    sequence: '/proj/phanstiel_lab/References/GENCODE.GRCh38.p10/Sequence/GRCh38.p10.genome.fa.gz' # Path to reference sequence gzipped fasta file. The associated sequence dictionary must be in the same directory as this file.

    ## RNA-seq read allele count thresholds 
    minTotalAlleleCounts: 10 # Minimum number of total read counts from both alleles of a variant to consider a donor as a heterozygote from RNA.
    minAlleleCounts: 2 # Minimum number of read counts from either allele of a variant to consider a donor as a heterozygote from RNA.
    ```

3. Submit the AI processing workflow to pre-process genotype data and RNA-seq data for AI analysis, combine allele counts, fill in missing data, and assign weights based on donor and genotyping concordance to all unique variants of the included samples.

    ```bash
    sbatch runAIprocessing
    ```
This workflow will produce the following key files:

- `output/AI/alleleCounts.csv`: Unfiltered variant allele counts with all merged and checked data for all donors and conditions. 
Columns of this file are `variantID`, `refCount`, `altCount`, `donor`, `condition`, and `weight`.
- `output/AI/numVariantHets.csv`: Quantification of number of donor heterozygotes for each variant.
- `output/AI/alleleCountsplits.txt`: A file listing the split files for alleleCounts.csv, which can be found in the `alleleSplits`
subfolder. These split files are for use in the final filtering workflow. 


4. Edit `config/config_filterVariants.yaml` for the third part of the allelic imbalance processing workflow with the following parameters:

    ```yaml
    ## Path to csv file that quantifies the number of heterozygote donors for each variant from AIprocessing workflow.
    numVariantHets: 'output/AI/numVariantHets.csv'

    ## Threshold of number of heterozygotes a variant for filtering.
    minHets: 5
    ```

5. Submit the third workflow to filter allele counts based on a heterozygote threshold and pivot data into an allele counts matrix and
weights matrix for use in `DESeq2`.

    ```bash
    sbatch runFilterVariants
    ```
This workflow will clean split files and produce the following key files:

- `output/AI/alleleCountsMatrix.csv`: The allele counts matrix to be used in `DESeq2`.
- `output/AI/weightsMatrix.csv`: The weights matrix to be used in `DESeq2`.
- `output/AI/colData.csv`: The data describing the columns (donor, condition, and allele) in the allele counts and weights
matrices.
