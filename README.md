# Easy Downstream Analysis of RNA Seq Data Package (easydars) 
Downstream analysis of bulk RNA sequencing data that has been run through DESeq2, including the creation of volcano plots, heat maps, PCA plots, and preparation of ranked genes lists and .gmt files for Fast Gene Set Enrichment Analysis.

## How to Install
```
install.packages("devtools")
library(devtools)
devtools::install_github("eitm-org/easydars")
library(easydars)
```

## Usage
### `DE_pca`
Inputs:
- DESeq dataset (dds)
- name of coefficient/design of experiment (xvar)
- number of top features to use (default 500)
- the PC to be on x-axis (default 1)
- the PC to be on y-axis (default 2)

Outputs:
A principal component analysis plot 

### `DE_volcano`
Inputs:
- counts data that contains both a gene_name and gene_id column (count_df)
- DESeq dataset (dds)
- DESeq results table (res)
- name of coefficient/design of experiment (xvar)
- name of your comparison 
- the names of the groups being compared
- shrinkage method (apeglm and ashr recommended(

Outputs:
A volcano plot
