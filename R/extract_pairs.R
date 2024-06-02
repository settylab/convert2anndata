#' Extract Pairwise Data from SingleCellExperiment
#'
#' This function extracts pairwise data (colPairs or rowPairs) from a SingleCellExperiment
#' object and returns it as a list of sparse matrices.
#'
#' @param pairs_function A function to extract pairwise data (e.g., colPairs or rowPairs).
#' @param sce A SingleCellExperiment object from which to extract the pairwise data.
#' @return A list of pairwise data as sparse matrices.
#' @importFrom SingleCellExperiment colPairs rowPairs
#' @importFrom Matrix sparseMatrix
#' @importFrom methods is
#' @export
extract_pairs <- function(pairs_function, sce) {
  pairs_list <- list()
  
  if (length(pairs_function(sce)) > 0) {
    pairs <- pairs_function(sce)
    for (name in names(pairs)) {
      if (is(pairs[[name]], "SelfHits")) {
        indices <- as.matrix(pairs[[name]])
        pairs_list[[name]] <- sparseMatrix(
          i = indices[, 1],
          j = indices[, 2],
          x = rep(1, nrow(indices)),
          dims = c(ncol(sce), ncol(sce))
        )
      } else {
        pairs_list[[name]] <- as(pairs[[name]], "dgCMatrix")
      }
    }
  }
  
  return(pairs_list)
}
