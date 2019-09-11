# Differential gene expression, GO enrichment and correlation analysis of a randomized double-blind placebo-controlled study with Enriched marine oils supplements 

# Overview: 

To establish the relationship between enriched marine oil supplement dose, peripheral specialized pro-resolving mediators (**SPM**) and cellular responses, a randomized double-blind placebo-controlled study was conducted. From blood, SPM (LC-MS/MS-based profiling) and gene expression (RNA Seq) profiles were obtained.  

This repository contains the **R scripts** used to run differential gene expression (DEG) and Gene Ontology (GO) enrichment analysis using the RNA Seq data. It also contains the script used to run **Spearman correlations** between SMP concentrations vs cell response.

**NOTE:** RNA seq raw file (fastq files) were uploaded to the [**Qiagen Data Analysis Center**](www.qiagen.com/GeneGlobe), and were processed with the **QIASeq UPX Primary Analysis Tool** for the trimming of the reads, demultiplexing of the samples, the alignment of the reads to the reference genome and the count of unique molecular identifiers (UMI). **Pathway analysis** were performed uploading the differential expressed genes obtained from the previous analysis in [**KEGG**](https://www.genome.jp/kegg/) and [**Reactome**](https://reactome.org/) pathway databases.  

# System Requirements: 

## Hardware requirements: 

All the scripts and software used in the paper *Enriched marine oils supplements increase peripheral blood SPM concentrations and reprogram host immune responses: A randomized double-blind placebo-controlled study* were run in a standard computer (RAM: 8GB, CP$: 4 cores, 3.60 GHZ/core) with a maximum runtime of approx. 8 minutes for the more demanding script [2_GO_analysis_(TopGO)](https://github.com/eagomezc/2019_DGE_Correlation_Oil_supplements_study/blob/master/b_R_Scripts/2_GO_analysis_(TopGO).R). 

A computer with lower specs (e.g. 2GB of RAM) will work but some scripts will take longer to run. 

## System requirements:

All the R scripts were created and used on **Windows 10**:

**R version**: 3.5.1 

**R Studio version**: 1.1.456

The scripts should be compatible with Mac and Linux operating systems. 

For installing R and R Studio, follows the installation instructions [here](https://www.stats.bris.ac.uk/R/) and [here](https://www.rstudio.com/products/rstudio/download/). With typical installation times on a normal computer not exceeding 3h.

## Required R packages (libraries): 

The required packages to run all the scripts contained in this repository can be installed as follow: 

```
# Package edgeR and TopGO from Bioconductor:
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("edgeR")
BiocManager::install("TopGO")

# Packages ggplot2, gplots, ggrepel, dvtools and corrplot:
install.packages(c('ggplot2', 'gplots', 'ggrepel', 'dvtools', 'corrplot'))
```
# Content:

This repository contains three folders: 

## [a_Data](https://github.com/eagomezc/2019_DGE_Correlation_Oil_supplements_study/tree/master/a_Data)

This folder contains, separated by subfolders, the different file formats that has to be used to run the different R scripts. Each subfolder has the name of the specific script where they can be used, in addition to the number of the script, to make more clear what file is used in what folder. At the moment to download this repository in a local computer, it's important to remember that all the **input pathways in the scripts has to be changed**.

The subfolders are:

**1_DGE_analysis_(Edge_R)**: Contains a tab-delimited table with the raw RNA-seq read counts (Patients as columns and genes as rows) and a tab-delimited table with the classification information (group, batch, subject, sample quality) of the patients. 

**2_GO_analysis_(TopGO)**: Contains a tab-delimited table with all the genes and their respective adjust *p* values from the differential gene expression analysis and a tab-delimited table with all the genes and their respective associated GO terms for a specific genome.

**3_Correlation_analysis_(Corrplot)**: Contains a tab-delimited table with a column with sample's name, a column with the condition, and then a series of columns with the different SPM concentrations and cell responses. 

More details about the format of this files can be seen in the comments of each script. 

## [b_R_Scripts](https://github.com/eagomezc/2019_DGE_Correlation_Oil_supplements_study/tree/master/b_R_Scripts)

This folder contains the scripts used to run the DGE, the GO enrichment analysis and the Spearman correlations made in this study. 

The scripts are: 

**1_DGE_analysis_(Edge_R).R**: Using RNA-seq raw read counts, this script performs differential gene expression analysis using the package **Edge R**, that uses the quasi-likelihood method to identify differences in the expression levels between the placebo group and the fish oil supplement group. The script generates three files: the counts per million of each gene, the DGE results and the resulting volcano plot from the DGE analysis. 

**2_GO_analysis_(TopGO).R**: Using the adjust *p* values from the DGE analysis and the associated GO terms for each gene, this script runs a GO enrichment analysis, which identifies the common GO terms associated with the genes that are up or downregulated in the supplement group when compared with the placebo group, using the package **TopGO**. The scripts generates two files for each GO category (Biological process, Cellular component and Molecular function): one file with the distribution plot of each GO term in the significant genes subset and a table cointaing the summary of the statistically significant GO terms. The GO enrichment analysis is made using the **Kolmogorov–Smirnov test** (KS). 

**3_Spearman_correlations_(Corrplot).R**: This script is designed to make a nonparametric correlation (Spearman correlation) between different features; specifically, the script correlates SPM concentrations with cell response and plots the results using the **Corrplot** package. The output of this script are the correlation matrix plots and tables with the correlation and adjust *p* values for all the correlations. 

More details of how the scripts works can be seen in the comments of each script. 

## [c_Expected_Output](https://github.com/eagomezc/2019_DGE_Correlation_Oil_supplements_study/tree/master/c_Expected_Outputs)

This folder contains, separated by subfolders, the different expected outputs that can be obtain after running the R scripts. Each subfolder has the name of the specific script that generates it, in addition to the number of the script, to make more clear what file is the result of what script. At the moment to download this repository in a local computer, it's important to remember that all the **output pathways in the scripts has to be changed**.

The subfolders are:

**1_DGE_analysis_(Edge_R)**: The expected results from this script are a tab-delimited file containing a table with the gene's names and the count per million (CPM) of each gene in the different samples (one sample per column); a tab-delimited table with all the the gene's names, their log(FC), log(CPM), F value, *p* value and adjust *p* value (FDR) from the DGE analysis, and the volcano plot from these results. 

**2_GO_analysis_(TopGO)**: The expected results from this script, for all the GO categories, is a tab-delimited file containing a table with the significant GO terms, the number of times it was annotated, the number of times it was annotated in significant genes from DGE, the expected number of times to be annonated, the *p* values and adjust *p* values (FDR) from the KS test and the significant genes associated with it. A second file, a pdf file, contains the distribution plots of the GO terms in all the genes and the significant genes; with this file is possible check which GO terms are actually showing a tendency of be annotated only in the significant genes from DGE. 

**3_Correlation_analysis_(Corrplot)**: The expected results from this script is a correlation matrix plot showing positive and negative correlation in different colors (blue for positive correlations and red for negative correlations). Besides that, the script produces the correlation values from the matrix (tab-delimited table) and their associated adjust *p* value (FDR) (tab-delimited table).

More details about how this files are generated can be seen in the comments of each script. 

# Reference: 

 If you use these codes for your research please cite these paper:
 
 *To be provided*
