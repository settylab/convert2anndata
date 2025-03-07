#' Process Main Assay
#'
#' This function processes the main assay in the SingleCellExperiment object to be used as the primary data matrix in the AnnData object.
#'
#' @param sce A SingleCellExperiment object.
#' @param assayName The name of the assay to use as the main data matrix in AnnData.
#' @return The main data matrix (anndata.X) for the AnnData object.
#' @importFrom SummarizedExperiment assay
#' @importFrom Matrix t
#' @export
process_main_assay <- function(sce, assayName) {
  timestamped_cat(sprintf(
    "Processing assay '%s' for anndata.X...\n",
    assayName
  ))
  alt_assayName <- assayName
  if (!(assayName %in% names(assays(sce)))) {
    alt_assayName <- names(assays(sce))[1]
    timestamped_cat(sprintf(
      "WARNING: The specified assay '%s' is not available. Using '%s' as active layer instead.\n",
      assayName, alt_assayName
    ))
  }
  X <- Matrix::t(assay(sce, alt_assayName))
  # Ensure X is in CsparseMatrix format
  X <- methods::as(X, "CsparseMatrix")
  timestamped_cat(sprintf(
    "Using '%s' assay as the main data matrix (converted to CsparseMatrix).\n",
    alt_assayName
  ))
  return(X)
}
