#' Extract Counts Matrix from Seurat Assay
#'
#' Retrieves the counts matrix from the default assay of the Seurat object.
#'
#' @param assay_object The default assay object from the Seurat object.
#' @return A counts matrix.
#' @importFrom Seurat GetAssayData
#' @export
extract_counts_matrix <- function(assay_object) {
  counts_matrix <- tryCatch(
    {
      if (inherits(assay_object, "Assay5")) {
        if ("counts" %in% Layers(assay_object)) {
          GetAssayData(assay_object, layer = "counts")
        } else {
          # If counts layer doesn't exist, try default
          GetAssayData(assay_object)
        }
      } else if (inherits(assay_object, "Assay")) {
        GetAssayData(assay_object, slot = "counts")
      } else {
        stop("Unsupported assay class: ", class(assay_object))
      }
    },
    error = function(e) {
      stop("Failed to retrieve counts matrix: ", e$message)
    }
  )

  if (nrow(counts_matrix) == 0 || ncol(counts_matrix) == 0) {
    stop("The counts matrix is empty. Ensure the Seurat object contains data.")
  }

  # Ensure that counts_matrix has rownames and colnames
  if (is.null(rownames(counts_matrix))) {
    rownames(counts_matrix) <- rownames(assay_object)
  }
  if (is.null(colnames(counts_matrix))) {
    colnames(counts_matrix) <- colnames(assay_object)
  }

  return(counts_matrix)
}
