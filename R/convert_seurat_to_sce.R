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
