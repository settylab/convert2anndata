#' Process Other Assays
#'
#' This function processes assays other than the main assay in the SingleCellExperiment object.
#'
#' @param sce A SingleCellExperiment object.
#' @param assayName The name of the assay used as the main data matrix.
#' @return A list of other assays for the AnnData object.
#' @importFrom SummarizedExperiment assays
#' @importFrom Matrix t
#' @export
process_other_assays <- function(sce, assayName) {
  timestamped_cat(sprintf("Processing assays other than '%s'...\n", assayName))
  all_assays <- assays(sce)
  all_assays <- all_assays[!names(all_assays) %in% assayName]
  all_assays <- lapply(all_assays, Matrix::t)
  timestamped_cat("Assays processed.\n")
  return(all_assays)
}
