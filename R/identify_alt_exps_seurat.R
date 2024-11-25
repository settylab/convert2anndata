#' Identify Assays to Move to altExps in Seurat Object
#'
#' Identifies assays with mismatching cells or features and prepares them to be moved to `altExps`.
#'
#' @param data The Seurat object.
#' @return A vector of assay names to be moved to `altExps`.
#' @importFrom Seurat GetAssayData
#' @export
identify_alt_exps_seurat <- function(data) {
  assays_to_include_in_altExp <- setdiff(names(data@assays), c("RNA", "RNA_combined"))

  # Filter out empty assays
  assays_with_data <- assays_to_include_in_altExp[vapply(assays_to_include_in_altExp, function(assay_name) {
    assay <- data@assays[[assay_name]]

    # Handle Assay5 separately from older Assay objects
    if (inherits(assay, "Assay5")) {
      counts <- GetAssayData(assay, layer = "counts")
      has_data <- (nrow(counts) > 0 && ncol(counts) > 0) || length(assay@layers) > 0
    } else if (inherits(assay, "Assay")) {
      counts <- GetAssayData(assay, slot = "counts")
      has_data <- nrow(counts) > 0 && ncol(counts) > 0
    } else {
      # Unsupported assay type
      has_data <- FALSE
    }

    return(has_data)
  }, FUN.VALUE = logical(1))]

  return(assays_with_data)
}
