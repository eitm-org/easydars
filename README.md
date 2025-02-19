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
- dds: DESeq dataset
- xvar: Coefficient or design of experiment (column in metadata)
- top_n: Number of top features to use (default 500)
- PC_of_interest_x: PC to be on x axis (default 1)
- PC_of_interest_y: PC to be on y axis (default 2)
  
Outputs:
- A principal component analysis plot 

### `DE_volcano`
Inputs:
- count_df: RNA Seq raw counts dataframe, including columns gene_id and gene_name
- dds: DESeq dataset
- res: DESeq2 results (after running results() on DESeq dataset)
- shrinkage_method: Shrinkage method (apeglm, ashr, or normal)
- xvar: Coefficient or design of experiment
- name_of_comparison: Name of comparison (column in metadata)
- first_group: Name of first group of comparison
- second_group: Name of second group of comparison
  
Outputs:
- A volcano plot

### `DE_heat_map`
Inputs:
- count_df: RNA Seq raw counts dataframe, including columns gene_id and gene_name
- col_data: Column data (sample metadata)
- dds: DESeq dataset
- res: DESeq2 results (after running results() on DESeq dataset)
- xvar: Coefficient or design of experiment
- name_of_comparison: Name of comparison
- top_n: Number of genes to plot (default 50)
- genes_of_interest: List of genes of interest (default none)
  
Outputs: 
- A heat map

### `prepare_gmt`
Inputs:
- gmt_file: gmt files containing genes sets
- genes_in_data: A list of genes present in the data

Outputs:
- A prepared .gmt file

### `prepare_rankings`
Inputs:
- log2foldchange: Log2 fold change of genes
- pvalues: P values of genes
- row_names: Gene names

Outputs:
- A ranked list of genes
