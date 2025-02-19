#' Creates a heat map
#' @param count_df RNA Seq raw counts dataframe, including columns gene_id and gene_name
#' @param col_data Column data (sample metadata)
#' @param dds DESeq dataset
#' @param res DESeq2 results (after running results() on DESeq dataset)
#' @param xvar Coefficient or design of experiment
#' @param name_of_comparison Name of comparison
#' @param top_n Number of genes to plot (default 50)
#' @param genes_of_interest List of genes of interest (default none)
#' @return A heat map
#' @export

DE_heat_map <- function(count_df,
                        col_data,
                        dds,
                        res,
                        xvar,
                        name_of_comparison,
                        top_n = 50,
                        genes_of_interest = NULL
                        ) {
  # Order results by adjusted p-value
  resOrdered <- res[order(res$padj),]

  # Determine genes to use for the heatmap
  if (!is.null(genes_of_interest)) {
    if(all(genes_of_interest %in% count_df$gene_name)){
      genes_to_plot <- count_df %>%
        filter(gene_name %in% genes_of_interest) %>%
               pull(gene_id)
      }
    else if(genes_of_interest %in% count_df$gene_id){
      genes_to_plot <- genes_of_interest
      }
    else{
      warning("genes_of_interest contains invalid gene names or ids.")
      }
  }
  else {
    resSig <- subset(resOrdered, padj < 0.05)  # Filter for significant genes
    genes_to_plot <- rownames(head(resSig, n = top_n))  # Select top_n genes by p-value
  }

  # Extract normalized counts for the selected genes
  normalized_counts <- vst(dds, blind = FALSE)  # Perform variance stabilizing transformation
  mat <- assay(normalized_counts)[genes_to_plot,]

  # Map Ensemble IDs to gene names using count_df
  res_df <- as.data.frame(resOrdered) %>%
    rownames_to_column(var = "gene_id") %>%
    left_join(count_df, by = "gene_id") %>%  # Join with count_df using gene_id
    mutate(gene_name = ifelse(is.na(gene_name), gene_id, gene_name)) %>%  # Replace NA with gene_id if no match found
    relocate(gene_name)  # Relocate gene_name to the front


  # Replace rownames with gene names
  rownames(mat) <- res_df$gene_name[match(rownames(mat), res_df$gene_id)]

  # Use sample group information for column names
  colnames(mat) <- col_data[[xvar]] # name_of_comparison is a string

  # Create the heatmap
  pheatmap(
    mat,
    cluster_rows = TRUE,    # Cluster genes (rows)
    cluster_cols = TRUE,    # Cluster samples (columns)
    show_rownames = TRUE,   # Show gene names
    show_colnames = TRUE,   # Show sample names
    fontsize_row = 8,       # Adjust row font size
    fontsize_col = 10,      # Adjust column font size
    legend = TRUE,          # Show color legend
    color = colorRampPalette(c("blue", "white", "red"))(500)  # Color scale (blue to red)
  )
}
