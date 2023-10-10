# Single-cell RNA Sequencing Analysis of Lung Cells in COVID-19 Patients with Diabetes, Hypertension, and Comorbid DiabetesHypertension

# Table of Contents
- [What is this?](#what-is-this)
- [How can I use this data, and where can I find it?](#how-can-i-use-this-data-and-where-can-i-find-it)
	- [Downloading Data files](##downloading-data-files)
		-[Single cell RNA sequencing datasets](###Single-cell-RNA-sequencing-datasets)
   		-[Bulk RNA sequencing datasets](####-Bulk-RNA-sequencing-datasets)
- [Analysis and visualization programs](#analysis-and-visualization-programs)
	- [R and R's integrated developmental environment RStudio:](#r-and-rs-integrated-developmental-environment-rstudio)
	- [scRNAseq analysis pipeline SEURAT developed by the Satija lab:](#scrnaseq-analysis-pipeline-seurat-developed-by-the-satija-lab)
- [cell–cell communication](##cell--cell-communication)
 	-[Perform NicheNet analysis starting from a Seurat object](###Perform-NicheNet-analysis-starting-from-a-Seurat-object)
  	-[CellChat](###CellChat)
- [Lead Contacts](#Lead-contacts)


## What is this?
This repository contains coding scripts utilized for the analysis performed in the "Single-cell RNA Sequencing Analysis of Lung Cells in COVID-19 Patients with Diabetes, Hypertension, and Comorbid DiabetesHypertension (Xin Zhang et. al, 2023). The purpose of providing the code here is to allow for transparency and robust data-analysis reproducibility. Most of the steps used for data analysis and visualization have been optimised for an average computing environment (for the year 2023). Some analyses however, require a high-performace computing environment (see computing environment). The methodology has already been described extensively in the manuscript. However, this analysis relies heavily on powerful scRNAseq analysis algorithms developed by the [Satija lab](https://satijalab.org/), namely [Seurat](https://satijalab.org/seurat/) [(Butler et al., 2018: Nature Biotechnology;](https://www.nature.com/articles/nbt.4096) [Stuart et al., 2018: Cell)](https://www.sciencedirect.com/science/article/pii/S0092867419305598?via%3Dihub) (for a complete list of dependencies and code utilized see analysis & visualization programs).

# How can I use this data, and where can I find it?
## Downloading Data files
#### Single cell RNA sequencing datasets
We downloaded the available COVID-19 and IPF-related scRNA-seq datasets from [Gene Expression Omnibus (GEO)](https://www.ncbi.nlm.nih.gov/geo/) and collected the datasets as follows:
1) COVID-19: The datasets included in our analysis were [GSE171524](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE171524), [GSE171668](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE171668), [GSE149878](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE149878), [GSE161382](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE161382), and [GSE163919](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE163919). 
2) Idiopathic pulmonary fibrosis (IPF): [GSE132771](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE132771), [GSE135893](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE135893), [GSE122960](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE122960), [GSE128033](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE128033), [GSE128169](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE128169), [GSE136831](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE136831), [GSE159354](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE159354).
#### Bulk RNA sequencing datasets
Data is part of the GSE47460 from GEO.
Contact lead author for Seurat object.

## Analysis and visualization programs
### R and R's integrated developmental environment RStudio
1、[R v4.2.2 (x64 bit)](https://cran.r-project.org/bin/windows/base/old/)

2、[Tutorial for R](https://cran.r-project.org/doc/manuals/r-release/R-intro.html)

3、[Tutorial for RStudio](https://resources.rstudio.com/) 


### scRNAseq analysis pipeline SEURAT developed by the Satija lab
[Source code for Seurat](https://cran.r-project.org/web/packages/Seurat/index.html)
[Tutorials for Seurat](https://satijalab.org/seurat/)

## cell-cell communication
### Perform NicheNet analysis starting from a Seurat object: step-by-step analysis（https://github.com/saeyslab/nichenetr/blob/master/vignettes/seurat_steps.md）
### CellChat
（https://github.com/sqjin/CellChat）
## Lead Contacts
Dr. Zhang,lanlan. Department of Respiratory and Critical Care Medicine, West China Hospital, Sichuan University. Chengdu, China. E-mail: lanlanzhangsc@gmail.com
