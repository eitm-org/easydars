
#' Creates a volcano plot
#' @param count_df RNA Seq raw counts dataframe, including columns gene_id and gene_name
#' @param dds DESeq dataset
#' @param res DESeq2 results (after running results() on DESeq dataset)
#' @param shrinkage_method Shrinkage method (apeglm, ashr, or normal)
#' @param xvar Coefficient or design of experiment
#' @param name_of_comparison Name of comparison (column in metadata)
#' @param first_group Name of first group of comparison
#' @param second_group Name of second group of comparison
#' @return A volcano plot
#' @export

DE_table_volcano <- function(count_df,
                             dds,
                             res,
                             shrinkage_method,
                             xvar,
                             name_of_comparison,
                             first_group,
                             second_group)
  {
  # Extract significant genes
  resOrdered <- res[order(res$padj),]
  resSig <- subset(resOrdered, padj < 0.05)
  resSig <- as.data.frame(resSig) %>%
    rownames_to_column(var = "gene_id") %>%
    left_join(count_df, by = "gene_id") %>%
    mutate(gene_name = ifelse(is.na(gene_name), gene_id, gene_name))
  resSig <- relocate(resSig, gene_name)
  min_pvalue <- min(resSig$pvalue[resSig$pvalue > 0], na.rm = TRUE)/10
  resSig$pvalue[resSig$pvalue == 0] <- min_pvalue
  resSig$padj[resSig$padj == 0] <- min_pvalue

  # Shrink the log fold changes if shrinkage_method is specified
  if (!missing(shrinkage_method)) {
    if (shrinkage_method == "apeglm") {
      res_shrink <- lfcShrink(dds,
                              coef = name_of_comparison,
                              res = res,
                              type = "apeglm")
    } else if (shrinkage_method == "ashr") {
      res_shrink <- lfcShrink(dds,
                              contrast = c(xvar, first_group, second_group),
                              res = res,
                              type = "ashr")
    }
    else if (shrinkage_method == "normal") {
      res_shrink <- lfcShrink(dds,
                              contrast = c(xvar, first_group, second_group),
                              res = res,
                              type = "normal")
    } else {
      stop("Invalid shrinkage method. Choose either 'apeglm', 'ashr', or 'normal'.")
    }

    # Convert to data frame and merge with gene names
    res_shrink_df <- as.data.frame(res_shrink) %>%
      rownames_to_column(var = "gene_id") %>%
      left_join(count_df, by = "gene_id") %>%
      mutate(gene_name = ifelse(is.na(gene_name), gene_id, gene_name))
    res_shrink_df <- relocate(res_shrink_df, gene_name)
    res_shrink_df$pvalue[res_shrink_df$pvalue == 0] <- min_pvalue
    res_shrink_df$padj[res_shrink_df$padj == 0] <- min_pvalue
  } else {
    res_shrink_df <- NULL
  }

  # Create the volcano plot
  volcano_plot <- EnhancedVolcano(res,
                                  lab = if (!is.null(res_shrink_df)) res_shrink_df$gene_name else resSig$gene_name,
                                  title = if (!is.null(first_group) && !is.null(second_group)) paste(first_group, "vs.", second_group) else "Volcano Plot",
                                  x = 'log2FoldChange',
                                  y = 'padj',
                                  xlab = bquote(~Log[2]~ 'fold change'),
                                  ylab = bquote(~-Log[10]~italic(padj)),
                                  legendLabels = c('Not significant',
                                                   bquote(~Log[2]~ 'fold change'),
                                                   'padj',
                                                   bquote('padj and' ~Log[2]~ 'fold change')))
  return(list(volcano_plot = volcano_plot,
              de_genes = resSig,
              res_shrink = res_shrink_df))
}



