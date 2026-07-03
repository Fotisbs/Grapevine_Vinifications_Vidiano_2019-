# ***Metataxonomic and metatranscriptomic analysis reveal microbial succession and metabolic pathways activated during spontaneous and inoculated vinification ***

### By Bekris F. <sup>1+</sup>, Papadopoulou E. <sup>1+</sup>, Vasileiadis S. <sup>1</sup>, Lola D. <sup>2</sup>, Kotseridis Y. <sup>2</sup>, Karpouzas D.G <sup>1*</sup>

### (\* corr. author)
### (\+ contributed equally to this work)

<sup>1</sup> University of Thessaly, Department of Biochemistry and Biotechnology, Laboratory of Plant and Environmental Biotechnology, 41500 Viopolis – Larissa, Greece

<sup>2</sup> Agricultural University of Athens, Department of Food Science and Human Nutrition, Laboratory of Enology and Alcoholic Drinks (LEAD),  11855 Iera Odos St. 75 – Athens, Greece


## Repository overview

This repository contains all scripts used to reproduce the microbiome and metatranscriptomic analyses presented in the manuscript. The repository includes data retrieval scripts, preprocessing workflows, statistical analyses and figure generation scripts.

To obtain the repository, install Git (if not already installed https://github.com/git-guides/install-git), open a terminal and clone the repository:

```
$ git clone https://github.com/Fotisbs/Grapevine_Vinifications_Vidiano_2019-.git
```

## Alternatively, the repository can be downloaded as a ZIP archive directly from GitHub.

Unless otherwise stated, all commands assume that the repository root directory ("Grapevine_Vinifications_Vidiano_2019-") is used as the working directory. The required sequencing datasets can be downloaded directly from the NCBI Sequence Read Archive using the scripts provided in each module.

## Repository structure
Fungi/
    0.DownloadData/
    1.Demultiplex/
    2.PhyloseqObjectPreparation/
    3.DataAnalysis/

Bacteria/
    0.DownloadData/
    1.Demultiplex/
    2.PhyloseqObjectPreparation/
    3.DataAnalysis/

Metatranscriptomic/
    0.DownloadData/
    1.FunctionalAnnotation/
    2.DataAnalysis/
	
## Description of the order of executed scripts.
# Microbiome (Metataxonomic) analyses

For Fungi and Bacteria files, steps 0-2 concern the data retrieval from NCBI and preprocessing (demultiplex) and phyloseq object construction, while step 3 and the subfolders concern the actual data analysis.

0) First, it is necessary to download the sequencing data.
To do so, you need to enter the "0.DownloadData" subfolder of "Fungi" and "Bacteria" folders accordingly and execute the "fetch_data.sh" bash script for batch (01), this assumes that you are located at the working directory "Grapevine_Vinifications_Vidiano_2019-".

The script downloads all raw amplicon sequencing reads deposited in the NCBI Sequence Read Archive using the corresponding SRR accession numbers listed in the batch files.Once the download is done, you need to combine all forward reads to a single file and all reverse reads to another file as well.
```
for i in {01}
do
	cd Fungi/0.DownloadData/batch${i}
	sh -x fetch_data.sh
	cat *_1.fastq | gzip > forward.fastq.gz
	cat *_2.fastq | gzip > reverse.fastq.gz
	cd ../../../
	cd Bacteria/0.DownloadData/batch${i}
	sh -x fetch_data.sh
	cat *_1.fastq | gzip > forward.fastq.gz
	cat *_2.fastq | gzip > reverse.fastq.gz
	cd ../../../
done
```

1) Then you need to demultiplex the data according to our own demultiplexing method using our in-house script.
This step requires Flexbar v3.0.3 together with the mapping file (map_file) provided in the corresponding folder. A detailed description of our in-house multiplexing approach is provided in our [previous work] (https://github.com/SotiriosVasileiadis/mconsort_tbz_degr#16s).
You need to enter the folder Fungi (or Bacteria)/1.Demultiplex and run the following commands (change the MY_PROCS variable to whatever number of logical processors you have available and want to devote),
the following commands are going to save the demultiplexed files in the Fungi(or Bacteria)/1.Demultiplex/demux_out folder.

## The demultiplexing workflow follows the protocol described in: https://github.com/SotiriosVasileiadis/mconsort_tbz_degr#16s

2) Following, the "Vinification Vidiano 2019 Quality-Classification-Phyloseq Object.R" script of the Fungi(or Bacteria)/2.PhyloseqObjectPreparation folder is run in order to prepare the final phyloseq object to be used in the data analysis described below. Before running the script make sure that the necessary reference databases are found in the same folder. The taxonomic annotations of the resulting fungal and bacterial ASVs were performed using the UNITE ITS v.8.2 (04.02.2020) (Morrison-Whittle et al., 2017) and the Silva v.138 (Yilmaz et al., 2014) databases as references respectively. The sample metadata file (samdf.txt), included in the repository, is also required for construction of the phyloseq objects.
```
cd Fungi/2.PhyloseqObjectPreparation
# fetch the databases
wget https://files.plutof.ut.ee/public/orig/1D/B9/1DB95C8AC0A80108BECAF1162D761A8D379AF43E2A4295A3EF353DD1632B645B.gz
# run the R script
Fungi Vinification Vidiano 2019 Quality-Classification-Phyloseq Object.r
cd ../../
cd Bacteria/2.PhyloseqObjectPreparation
# fetch the databases
wget https://zenodo.org/record/4587955/files/silva_nr99_v138.1_train_set.fa.gz
wget https://zenodo.org/record/4587955/files/silva_nr99_v138.1_wSpecies_train_set.fa.gz
tar vxf *.gz
# run the R script
Bacteria Vinification Vidiano 2019 Quality-Classification-Phyloseq Object.r
cd ../../
```
3) The DataAnalysis folder contains independent R scripts reproducing all microbiome analyses and figures presented in the manuscript.

```
3a. Taxonomic composition (bar plots)

3b. Beta-diversity (NMDS)

3c. Beta-diversity statistics (PERMANOVA)

3d. Rarefaction curves

3e. Alpha-diversity metrics

3f. Differential abundance heatmaps
```

# Metatranscriptomic analyses

Metatranscriptomic analyses are organized into three sequential steps: (0) retrieval of raw sequencing data, (1) functional annotation using the SAMSA2 pipeline, and (2) statistical analyses in R.

0) First, it is necessary to download the RNA sequencing data.
Download the RNA sequencing data from the NCBI Sequence Read Archive.
To do so, you need to enter the "0.DownloadData" subfolder of "Metatranscriptomic" and execute the "fetch_data.sh" bash script for batch (01), this assumes that you are located at the working directory "Grapevine_Vinifications_Vidiano_2019-". The NCBI submitted RNA sequences are includes at those batch/files.The script is based on the SRR accession numbers for each batch file and can be found in the 0.DownloadData folder as a.txt file.
Once the download is done, you need to combine all forward reads to a single file and all reverse reads to another file as well.
```
for i in {01}
do
    cd Metatranscriptomic/0.DownloadData/batch${i}
    sh -x fetch_data.sh
    cat *_1.fastq | gzip > forward.fastq.gz
    cat *_2.fastq | gzip > reverse.fastq.gz
    cd ../../../
done
```

1) Functional annotation with SAMSA2
Raw paired-end RNA sequencing reads were processed using the SAMSA2 v2.2.0 pipeline (Westreich et al., 2018).
The complete preprocessing workflow is provided in

```
master_script_final.bash
```

## The pipeline performs the following steps automatically:

# Quality trimming using Trimmomatic
# Paired-end read merging using PEAR
# Removal of residual rRNA sequences using SortMeRNA
# Functional annotation against the RefSeq database using DIAMOND BLASTx
# Functional annotation against the SEED Subsystems database
# Generation of organism-level, gene-level and subsystem count tables
# Initial differential abundance analysis using the SAMSA2 R scripts

The pipeline can be executed with

```
sh master_script_final.bash
```

2) Statistical analyses in R

The downstream statistical analyses presented in the manuscript were performed in R using the annotation tables generated by SAMSA2.

The required input files are

```
FunctionalAnnotation.txt
ExperimentalDesign.txt
```

The processed functional annotation table (FunctionalAnnotation.txt) and the experimental design file (ExperimentalDesign.txt) are provided in this repository and serve as the input for all downstream statistical analyses. The original SAMSA2 pipeline (master_script_final.bash) used to generate these annotation tables is also included for reproducibility.

The analysis folder contains separate scripts reproducing all transcriptomic figures presented in the manuscript:

2a. NMDS ordination

2b. PERMANOVA

2c. Differential gene expression (Volcano plot)

2d. Differential gene expression heatmaps

2e. Functional subsystem analysis

All downstream analyses were performed in R using DESeq2. Different statistical models were used depending on the analysis:

Volcano plot: differential expression was estimated using the model
```
design = ~ Vinification + Stage
```
to evaluate the overall effect of fermentation strategy while accounting for differences among fermentation stages.

Heatmaps of differentially expressed genes: genes were identified using pairwise contrasts between inoculated and spontaneous fermentations performed independently within each fermentation stage (Start, Middle and End). DESeq2-normalized counts were subsequently used for visualization.

Functional subsystem analysis: subsystem abundances were compared between inoculated and spontaneous fermentations independently within each fermentation stage using DESeq2 pairwise contrasts.


## Code Usage disclaimer<a name="disclaimer"></a>

The following is the disclaimer that applies to all scripts, functions, one-liners, etc. This disclaimer supersedes any disclaimer included in any script, function, one-liner, etc.

You running this script/function means you will not blame the author(s) if this breaks your stuff. This script/function is provided **AS IS** without warranty of any kind. Author(s) disclaim all implied warranties including, without limitation, any implied warranties of merchantability or of fitness for a particular purpose. The entire risk arising out of the use or performance of the sample scripts and documentation remains with you. In no event shall author(s) be held liable for any damages whatsoever (including, without limitation, damages for loss of business profits, business interruption, loss of business information, or other pecuniary loss) arising out of the use of or inability to use the script or documentation. Neither this script/function, nor any part of it other than those parts that are explicitly copied from others, may be republished without author(s) express written permission. Author(s) retain the right to alter this disclaimer at any time. This disclaimer was copied from a version of the disclaimer published by other authors in https://ucunleashed.com/code-disclaimer and may be amended as needed in the future.
