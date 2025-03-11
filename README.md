# :dna: Easy Downstream Analysis of RNA Seq Data Package (easydars) 
Downstream analysis of bulk RNA sequencing data that has been run through DESeq2, including the creation of volcano plots, heat maps, PCA plots, and preparation of ranked genes lists and .gmt files for Fast Gene Set Enrichment Analysis.

## How to Install
```
install.packages("devtools")
library(devtools)
devtools::install_github("eitm-org/easydars")
library(easydars)
```
## Usage
A full length tutorial can be found on this repo, titled easydars_tutorial.Rmd'

# DE_PCA
```
# returns Principal Component Analysis Plot from DESeq dataset and experimental design coefficient
DE_pca(dds = genotype_DEs, xvar = "sample_genotype")
```
![PCA Plot](/sample_plots/easydars_PCA.png "PCA")
# DE_table_volcano
```
# returns a volcano plot 
DE_table_volcano(count_df = count_df,
                 dds = genotype_DEs,
                 res = res,
                 shrinkage_method = "apeglm",
                 xvar = "sample_genotype",
                 name_of_comparison = "sample_genotype_WT_vs_KO",
                 first_group = "WT",
                 second_group = "KO"
                                                                )
```
![Volcano Plot](/sample_plots/easydars_VolcanoPlot.png "Volcano Plot")
# DE_heat_map
```
# returns a heat map of top 10 DE genes
DE_heat_map(count_df = count_df,
            col_data = col_data,
            dds = genotype_DEs,
            res = res,
            xvar = "sample_genotype",
            name_of_comparison = "sample_genotype_WT_vs_KO",
            top_n = 10
                                                             )

# returns a heat map of genes of interest
genes_of_interest = c("Prex1","Lrp3","Tspan6", "Mical1")
DE_heat_map(count_df = count_df,
            col_data = col_data,
            dds = genotype_DEs,
            res = res,
            xvar = "sample_genotype",
            name_of_comparison = "sample_genotype_WT_vs_KO",
            genes_of_interest = genes_of_interest
                                                             )
```
![Heatmap](/sample_plots/easydars_top10HeatMap.png "Heatmap")
![Heatmap](/sample_plots/easy_dars_genesHeatMap.png "Heatmap")
# prepare_rankings
```
# prepares a ranked list of genes for fgsea
prepare_rankings(log2foldchange = fc, pvalues = p_val, row_names = rownames)
```
# prepare_gmt
```
# prepares gmt files to be used in fgsea
pathways <- prepare_gmt(gmt_file = gmt_file, genes_in_data = my_genes)
```
