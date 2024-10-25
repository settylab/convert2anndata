#' Convert SingleCellExperiment to AnnData
#'
#' This function converts a SingleCellExperiment (SCE) object to an AnnData object.
#' It processes assays, dimensional reductions, and metadata, including alternative experiments (altExps).
#' The main assay is used as the primary data matrix in the AnnData object.
#'
#' @param sce A SingleCellExperiment object to be converted.
#' @param assayName The name of the assay to use as the main data matrix in AnnData. Defaults to "counts".
#' @param useAltExp Logical indicating whether to process and include alternative experiments (altExps). Defaults to TRUE.
#' @return An AnnData object containing the data from the SingleCellExperiment object.
#' @details The function first prints a summary of the input SingleCellExperiment object and checks for alternative experiments (altExps).
#' If altExps are found, they are processed recursively and stored in the `uns` slot of the AnnData object, unless the `useAltExp` parameter is set to FALSE.
#' The specified assay is used as the primary data matrix (anndata.X), and all other assays are stored in the `layers` slot.
#' Dimensional reductions are stored in the `obsm` slot, and metadata is organized into `obs`, `var`, `obsp`, `varp`, and `uns` slots as appropriate.
#' The function includes detailed error handling to ensure that all data is correctly transferred to the AnnData object.
#' @examples
#' \dontrun{
#' library(SingleCellExperiment)
#' library(anndata)
#' sce <- SingleCellExperiment::SingleCellExperiment(list(counts = matrix(1:4, ncol = 2)))
#' ad <- convert_to_anndata(sce)
#' }
#' @importFrom anndata AnnData
#' @export
convert_to_anndata <- function(sce, assayName = "counts", useAltExp = TRUE) {
  # Print a summary of the input SCE object
  timestamped_cat("Summary of SingleCellExperiment object:\n\n")
  print(sce)
  cat("\n")

  # Process altExperiments
  alt_exps_data <- process_alt_experiments(sce, assayName, useAltExp)
  sce <- alt_exps_data$sce
  alt_exps <- alt_exps_data$alt_exps

  # Process the main assay
  X <- process_main_assay(sce, assayName)

  # Process other assays
  layers <- process_other_assays(sce)

  # Process dimensional reductions
  obsm <- process_dimensional_reductions(sce)

  # Process obs and var data
  obs_data <- extract_data(colData, "obs/colData", sce)
  var_data <- extract_data(rowData, "var/rowData", sce)

  # Filter out reducedDims from obs_data
  available_reductions <- names(reducedDims(sce))
  reduction_prefixes <- paste0(tolower(available_reductions), "\\.")
  reduction_columns <- grep("^reducedDims\\.", names(obs_data$data), value = TRUE)
  if (length(reduction_columns) > 0) {
    obs_data$data <- obs_data$data[, !(names(obs_data$data) %in% reduction_columns), drop = FALSE]
  }

  # Ensure obs and var are correctly provided
  if (nrow(obs_data$data) != nrow(X) || ncol(obs_data$data) == 0) {
    obs_data$data <- NULL
  }

  if (nrow(var_data$data) != ncol(X) || ncol(var_data$data) == 0) {
    var_data$data <- NULL
  }

  # Retrieve colPairs and rowPairs data
  col_pairs_data <- extract_pairs(colPairs, "Column Pairs", sce)
  row_pairs_data <- extract_pairs(rowPairs, "Row Pairs", sce)

  # Report extracted pairwise data
  if (length(col_pairs_data) > 0) {
    timestamped_cat("Extracted colPairs data:\n", paste(names(col_pairs_data), collapse = ", "), "\n")
  } else {
    timestamped_cat("No colPairs data extracted.\n")
  }

  if (length(row_pairs_data) > 0) {
    timestamped_cat("Extracted rowPairs data:\n", paste(names(row_pairs_data), collapse = ", "), "\n")
  } else {
    timestamped_cat("No rowPairs data extracted.\n")
  }

  # Process metadata and pairwise matrices
  uns_data <- process_metadata_and_pairwise(sce, alt_exps, X)

  # Create a list of arguments for AnnData
  anndata_args <- list(
    X = X,
    layers = layers,
    obsm = obsm,
    varm = uns_data$varm,
    obsp = c(uns_data$obsp, col_pairs_data),
    varp = c(uns_data$varp, row_pairs_data),
    uns = uns_data$uns
  )

  if (!is.null(obs_data$data)) {
    anndata_args$obs <- obs_data$data
  }
  if (!is.null(var_data$data)) {
    anndata_args$var <- var_data$data
  }

  # Create AnnData object using the arguments list
  ad <- do.call(AnnData, anndata_args)

  timestamped_cat("Summary of the AnnData object:\n\n")
  print(ad)
  cat("\n")

  return(ad)
}
