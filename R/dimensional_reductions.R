#' Process Dimensional Reductions
#'
#' This function processes the dimensional reductions in the SingleCellExperiment object.
#'
#' @param sce A SingleCellExperiment object.
#' @return A list of dimensional reductions for the AnnData object.
#' @importFrom SingleCellExperiment reducedDims reducedDim
#' @export
process_dimensional_reductions <- function(sce) {
  timestamped_cat("Processing dimensional reductions...\n")
  available_reductions <- names(reducedDims(sce))
  if (length(available_reductions) > 0) {
    obsm <- lapply(available_reductions, function(rd_name) {
      reducedDim(sce, rd_name)
    })
    names(obsm) <- paste("X", tolower(available_reductions), sep = "_")
    timestamped_cat("Dimensional reductions processed.\n")
  } else {
    obsm <- list()
    timestamped_cat("No dimensional reductions found.\n")
  }
  return(obsm)
}
