#' Process Other Assays
#'
#' This function processes assays other than the main assay in the SingleCellExperiment object.
#'
#' @param sce A SingleCellExperiment object.
#' @param assayName (Optional) The name of the assay used as the main data matrix.
#' If not provided, all assays will be processed.
#' @return A list of processed assays.
#' @importFrom SummarizedExperiment assays
#' @importFrom Matrix t
#' @export
process_other_assays <- function(sce, assayName = NULL) {
  if (is.null(assayName)) {
    timestamped_cat("Processing all assays...\n")
    all_assays <- assays(sce)
  } else {
    timestamped_cat(sprintf("Processing assays other than '%s'...\n", assayName))
    all_assays <- assays(sce)
    all_assays <- all_assays[!names(all_assays) %in% assayName]
  }

  processed_assays <- lapply(all_assays, Matrix::t)
  timestamped_cat("Assays processed.\n")
  return(processed_assays)
}
