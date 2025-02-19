#' Prepares a ranked list of genes using a signed p-value metric to be used in fgsea
#' @param log2foldchange Log2 fold change of genes
#' @param pvalues P values of genes
#' @param row_names Gene names
#' @return A ranked list of genes
#' @export
prepare_rankings <- function(log2foldchange, pvalues, row_names){
  rankings <- sign(log2foldchange)*(-log10(pvalues)) # signed p value metric for ranking
  names(rankings) <- toupper(row_names) # capitalize gene names
  rankings <- rankings[!(names(rankings) == "")] %>% na.omit() # remove rankings without gene names
  rankings <- rankings[!duplicated(names(rankings))] # ensure no duplicates
  rankings <- sort(rankings,decreasing = TRUE) # sort rankngs from highest to lowest

  return(rankings)
}

#' Prepares gmt files to be used in fgsea
#' @param gmt_file gmt files containing genes sets
#' @param genes_in_data A list of genes present in the data
#' @return A prepared gmt file
#' @export


# prepare gmt files
prepare_gmt <- function(gmt_file, genes_in_data, savefile = FALSE) {
  # helper function - adjacency matrix to list
  matrix_to_list <- function(pws){
    pws.l <- list()
    for (pw in colnames(pws)) {
      pws.l[[pw]] <- rownames(pws)[as.logical(pws[, pw])]
    }
    return(pws.l)
  }
  my_genes <- toupper(genes_in_data)
  gmt <- gmtPathways(gmt_file)  # returns list of pathways from gmt file
  hidden <- toupper(unique(unlist(gmt)))
  # Convert gmt file to a matrix with genes as rows, GO annotations (columns), values as 0 or 1
  mat <- matrix(0, dimnames = list(hidden, names(gmt)),
                nrow = length(hidden), ncol = length(gmt))
  for (i in seq_along(gmt)) {
    mat[hidden %in% gmt[[i]], i] <- 1  # Set 1 if gene is in pathway
  }
  # Subset matrix to only include genes present in `genes_in_data`
  mat <- mat[intersect(rownames(mat), my_genes), ]
  # Filter gene sets (columns) that have more than 5 annotated genes
  mat <- mat[, colSums(mat) > 5, drop = FALSE]
  # Convert matrix to a list format if needed
  final_list <- matrix_to_list(as.data.frame(mat))
  if (savefile) {
    saveRDS(final_list, file = paste0(gsub('.gmt$', '', gmt_file), '_subset_', format(Sys.time(), '%d%m'), '.RData'))
  }
  return(final_list)
}
