#' Align Metadata to Cell Names in Seurat Object
#'
#' Ensures that the metadata is correctly aligned with the cell names in the counts matrix.
#'
#' @param data The Seurat object.
#' @param combined_counts The combined counts matrix.
#' @return A metadata DataFrame aligned with the counts matrix.
#' @importFrom Seurat Cells
#' @export
align_metadata_seurat <- function(data, combined_counts) {
  metadata <- tryCatch(
    {
      aligned_meta <- data@meta.data[colnames(combined_counts), , drop = FALSE]
      if (nrow(aligned_meta) != ncol(combined_counts)) {
        stop("Mismatch between cell names in counts matrix and colData.")
      }
      if (!all(rownames(aligned_meta) == colnames(combined_counts))) {
        stop("Mismatch between cell names in counts matrix and colData.")
      }
      aligned_meta
    },
    error = function(e) {
      stop("Error aligning metadata to cell names: ", e$message)
    }
  )
  return(metadata)
}
