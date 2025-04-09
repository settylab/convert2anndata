#' Ensure Sparse Matrix is in a Format Compatible with AnnData
#'
#' This function checks if a sparse matrix is in a format compatible with AnnData 
#' (CsparseMatrix format, e.g., dgCMatrix) and converts it if necessary.
#'
#' @param mat A matrix or Matrix object
#' @return The original matrix if already in a compatible format, or a converted
#'         matrix if conversion was needed
#' @importFrom methods is as
#' @importFrom Matrix Matrix
#' @export
ensure_csparse_matrix <- function(mat) {
  # Only process Matrix objects
  if (!methods::is(mat, "Matrix")) {
    return(mat)
  }
  
  # Check if already in CsparseMatrix format
  if (methods::is(mat, "CsparseMatrix")) {
    return(mat)  # Already in appropriate format (e.g., dgCMatrix)
  }
  
  # Convert to CsparseMatrix format if it's a different sparse matrix type
  if (methods::is(mat, "sparseMatrix")) {
    timestamped_cat("Converting matrix from ", class(mat)[1], " to CsparseMatrix format for AnnData compatibility\n")
    return(methods::as(mat, "CsparseMatrix"))
  }
  
  # Return original matrix for dense matrices
  return(mat)
}