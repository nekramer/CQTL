# Genopipe

A pipeline for QC'ing and preparing genotyping data for imputation, combining genotyping data with a reference panel to 
determine ancestry with EIGENSTRAT, and QC'ing genotype vcf files post-imputation.

## Converting files from GenomeStudio into PLINK format for use in this pipeline

### *Installing and Loading*

GenomeStudio is an Illumina software that is only made for Windows. If using a Mac, install a Windows mirror with the following steps:

1. Download the [Windows 10 Disc Image (ISO File)](https://www.microsoft.com/en-au/software-download/windows10ISO)
2. Open Applications > Utilities > **Boot Camp Assistant**
3. Follow the onscreen instructions to set up Boot Camp Assistant. When prompted, add downloaded Windows 10 ISO image.
4. Switch between Mac and Windows by restarting and holding **option**.

Installing GenomeStudio:

1. Download both **GenomeStudio** and **PLINK Input Report Plug-in** from Illumina [here](https://emea.support.illumina.com/array/array_software/genomestudio/downloads.html).

2. Follow onscreen instructions to download and setup GenomeStudio.

### *Loading data into GenomeStudio and exporting into PLINK format*
**Do not open the included sample sheet before loading into Genome Studio because it will mess up the barcode formatting.**

1. Start a new **Genotyping Project**. Select "Use sample sheet to load sample intensities."
2. On the next page of the Project Wizard, load the given sample sheet (.csv) under "Sample Sheet." Load the given **images directory** under "Data Repository." Load the given **Manifest directory** under "Manifest Repository."
3. On the final page of the Project Wizard, click the box to "Import cluster positions from a cluster file." Load the given **/Manifest > .egt file** as the cluster file. Click Finish.
4. To export a PLINK report, click **Analysis > Reports > Custom Report > PLINK**. A **map** and **ped** file should be generated. 

## Workflow
1. Clone repo into working directory:

    ```bash
    git clone --no-checkout https://github.com/nekramer/CQTL.git Genopipe
    cd Genopipe
    git sparse-checkout init --cone
    git sparse-checkout set Genopipe
    git checkout
    ```

2. Edit the comma-separated `geno.csv` with the batch name, PLINK map/ped file prefix, and the path to these files under the `Genotyping_Directory`
column. The names in the `Batch` columns will be used to create a group name for output files.

3. Edit `config/config.yaml` for parameters specific to your analysis. The parameters are described below:

    ```yaml
    ## Path to genotyping data samplesheet
    geno: 'geno.csv'

    ## Thresholds for variant inclusion in PLINK
    miss_call: 0.1 # PLINK --geno: Missing rate per SNP; A value of 0.1 will include only SNPs with a 90% genotyping rate (10% missing).
    maf: 0.01 # PLINK --maf: Minor allele frequency; A value of 0.01 will only include SNPs with a minor allele frequencey >= 0.01.
    hwe: .000001 # PLINK --hwe: Hardy-Weinberg equilibrium; Exclude markers that fail the Hardy-Weinberg test at a specified significance.
    remove: 'remove.txt' # Text file listing which, if any, samples to remove from data and which 
    batch they should be removed from.

    ## Genome-specific reference parameters
    ref: '/proj/phanstiel_lab/References/genomes/1000G/GRCh37/1000G_phase3_chrALL_biallelic' # Path to population reference data in PLINK binary format (.bed, .bim,       .fam files), all autosomes merged.
    panel: '/proj/phanstiel_lab/References/genomes/1000G/GRCh37/1000G_phase3.panel' # Path to panel file of above population reference data. Must have the columns         `sample`, `pop`, `super_pop`, and `gender`.
    sequence: '/proj/phanstiel_lab/References/genomes/hg19/Sequence/hg19.fa' # Path to reference sequence fasta file.

    ## Additional options
    pop_name: 'CQTL' # Name of data 'population' to use as the label in ancestry PCA plot.
    ```
To check the genome build of raw genotyping data obtained via a GenomeStudio project, check the `GenomeBuild` column of the .csv file
found under the Manifest directory.

4. Submit the first workflow with `sbatch`:

    ```bash
    sbatch runGenopipe
    ```

After running these steps the pipeline will produce the following key files:
- `output/ancestry/ancestry.pdf`: A plot of PC1 vs. PC2 of data merged with reference, colored by population.
- `output/imputation/{group}_chr{chr}.recode.vcf.gz`: Gzipped vcf files, separated by chromosome, which are prepared for impututation.

5. Impute data with imputation method of choice (i.e. Michigan Imputation Server) using an appropriate reference file.

6. Edit the comma-separated `vcf.csv` with the chromosome and the full path/filename to the imputed vcf of that chromosome in the `vcf_path` field.

7. Submit the second workflow to combine and QC these vcfs with `sbatch`:
    ```bash
    sbatch runGenopipeImputation
    ```
After running this workflow the pipeline will produce one qc'd, gzipped vcf file:
- `output/vcf/{group}_ALL_qc.vcf.gz`


