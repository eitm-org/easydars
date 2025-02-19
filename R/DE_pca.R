#' Creates a Principal Component Analysis (PCA) plot
#' @param dds DESeq dataset
#' @param xvar Coefficient or design of experiment (column in metadata)
#' @param top_n Number of top features to use (default 500)
#' @param PC_of_interest_x PC to be on x axis (default 1)
#' @param PC_of_interest_y PC to be on y axis (detaulf 2)
#' @return A PCA
#' @export

DE_pca <- function(dds,
                   xvar,
                   top_n = 500,
                   PC_of_interest_x = 1,
                   PC_of_interest_y = 2) {

  # Apply variance stabilizing transformation (VST) to the data
  vst_data <- vst(dds, blind = TRUE)

  # Extract the normalized counts
  counts_matrix <- assay(vst_data)

  # Get the top n features by variance
  var_genes <- order(apply(counts_matrix, 1, var), decreasing = TRUE)[1:top_n]
  top_counts <- counts_matrix[var_genes, ]

  # Perform PCA on the top n features
  pca_res <- prcomp(t(top_counts))

  # Extract the PCA data
  pca_data <- as.data.frame(pca_res$x)
  pca_data$xvar <- colData(dds)[[xvar]]

  # Plot the selected PCs (PC1 vs PC2 by default)
  ggplot(pca_data, aes_string(x = paste0("PC", PC_of_interest_x),
                              y = paste0("PC", PC_of_interest_y),
                              color = "xvar")) +
    geom_point(size = 3) +
    xlab(paste("PC", PC_of_interest_x)) +
    ylab(paste("PC", PC_of_interest_y)) +
    theme_minimal() +
    scale_color_discrete(name = xvar)
}
