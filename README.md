# Differential gene expression, GO enrichment and correlation analysis of a randomized double-blind placebo-controlled study with Enriched marine oils supplements 

# Overview: 

To establish the relationship between enriched marine oil supplement dose, peripheral specialized pro-resolving mediators (**SPM**) and cellular responses, a randomized double-blind placebo-controlled study was conducted. From blood, SPM (LC-MS/MS-based profiling) and gene expression (RNA Seq) profiles were obtained.  

This repository contains the **R scripts** used to run differential gene expression (DEG) and Gene Ontology (GO) enrichment analysis using the RNA Seq data. It also contains the script used to run **Spearman correlations** between SMP concentrations vs cell response, and SMP concentrations and marine oil supplement dose. 

**NOTE:** RNA seq raw file (fastq files) were uploaded to the [**Qiagen Data Analysis Center**](www.qiagen.com/GeneGlobe), and were processed with the **QIASeq UPX Primary Analysis Tool** for the trimming of the reads, demultiplexing of the samples, the alignment of the reads to the reference genome and the count of unique molecular identifiers (UMI). **Pathway analysis** were performed uploading the differential expressed genes obtained from the previous analysis in [**KEGG**](https://www.genome.jp/kegg/) and [**Reactome**](https://reactome.org/) pathway databases.  

# System Requirements: 

## Hardware requirements: 

All the scripts and software used in the paper *Enriched marine oils supplements increase peripheral blood SPM concentrations and reprogram host immune responses: A randomized double-blind placebo-controlled study* were run in a standard computer (RAM: 8GB, CP$: 4 cores, 3.60 GHZ/core) with a maximum runtime of approx. 3 minutes for the more demanding script [2_GO_analysis_(TopGO)](). 

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

# Packages ggplot2, gplots, ggrepel, dvtools, ggpubr and corrplot:
install.packages(c('ggplot2', 'gplots', 'ggrepel', 'dvtools', 'ggpubr', 'corrplot'))
```
# Content:



