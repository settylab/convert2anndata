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
  layers <- process_other_assays(sce, assayName)
  
  # Process dimensional reductions
  obsm <- process_dimensional_reductions(sce)
  
  # Process obs and var data
  obs_data <- extract_data(colData, "obs/colData", sce)
  var_data <- extract_data(rowData, "var/rowData", sce)
  
  # Ensure obs and var are correctly provided
  if (nrow(obs_data$data) != nrow(X) || ncol(obs_data$data) == 0) {
    obs_data$data <- NULL
  }
  
  if (nrow(var_data$data) != ncol(X) || ncol(var_data$data) == 0) {
    var_data$data <- NULL
  }

  # Process metadata and pairwise matrices
  uns_data <- process_metadata_and_pairwise(sce, alt_exps, X)

  # Create a list of arguments for AnnData
  anndata_args <- list(
    X = X,
    layers = layers,
    obsm = obsm,
    varm = uns_data$varm,
    obsp = uns_data$obsp,
    varp = uns_data$varp,
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

#' Convert Seurat or Other Object to SingleCellExperiment
#'
#' This function determines the class of a loaded object and converts it to a SingleCellExperiment object if necessary.
#' It handles Seurat objects, updating old Seurat v2 objects if detected, and converts them to SingleCellExperiment.
#' If the input object is already a SingleCellExperiment, it is returned as is.
#' For other object types, an attempt is made to convert them to SingleCellExperiment.
#'
#' @param data The loaded object to be converted.
#' @return A SingleCellExperiment object.
#' @details The function first checks if the input object is a Seurat object. If it is, it prints a summary of the object
#' and checks for indicators of Seurat v2, updating the object if necessary. It then converts the Seurat object to a SingleCellExperiment.
#' If the object is already a SingleCellExperiment, it is returned directly. For other object types, an attempt is made to convert them to
#' SingleCellExperiment.
#' @examples
#' \dontrun{
#' library(Seurat)
#' library(SingleCellExperiment)
#' seurat_obj <- CreateSeuratObject(counts = matrix(1:4, ncol = 2))
#' sce <- convert_seurat_to_sce(seurat_obj)
#' }
#' @export
convert_seurat_to_sce <- function(data) {
  object_class <- class(data)
  timestamped_cat("Input object class:", paste(object_class, collapse = ", "), "\n")

  if ("seurat" %in% tolower(object_class)) {
    timestamped_cat("Summary of input Seurat object:\n\n")
    suppressPackageStartupMessages(print(data))
    cat("\n")

    # Use tryCatch to safely check for Seurat v2 or v3 indicators
    raw_data <- tryCatch({
      if (!is.null(data@raw.data)) data@raw.data else NULL
    }, error = function(e) NULL)

    if (!is.null(raw_data)) {
      timestamped_cat("Old Seurat v2 object detected, attempting to update...\n")
      data <- Seurat::UpdateSeuratObject(data)
    }

    # Convert to SingleCellExperiment
    sce <- tryCatch({
      Seurat::as.SingleCellExperiment(data)
    }, error = function(e) {
      # Handle the case where layers might be empty
      Seurat::as.SingleCellExperiment(data, assay = NULL)
    })
  } else if ("SingleCellExperiment" %in% object_class) {
    sce <- data
  } else {
    suppressPackageStartupMessages({
      sce <- as(data, "SingleCellExperiment")
    })
  }
  return(sce)
}
