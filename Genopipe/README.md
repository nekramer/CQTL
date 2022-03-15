# Genopipe

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