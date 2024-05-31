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
  timestamped_cat(sprintf("Processing assay '%s' for anndata.X...\n", 
                          assayName))
  if (!(assayName %in% names(assays(sce)))) {
    timestamped_cat(
      sprintf(
        "Error: The specified assay '%s' is not available in the provided ",
        "SingleCellExperiment object.",
        assayName
      ),
      "Use -a or --assay to specify an available assay.\n"
    )
    timestamped_cat("Available assays are:", 
                    paste(names(assays(sce)), collapse = ", "), "\n")
    quit(status = 1, save = "no")
  }
  X <- Matrix::t(assay(sce, assayName))
  timestamped_cat(sprintf("Using '%s' assay as the main data matrix.\n", 
                          assayName))
  return(X)
}
