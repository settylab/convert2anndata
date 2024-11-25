#' Create Combined Assay in Seurat Object
#'
#' Creates a new assay in the Seurat object using the combined counts matrix.
#'
#' @param data The Seurat object.
#' @param combined_counts The combined counts matrix.
#' @param condition_labels Labels indicating the origin of each cell.
#' @return The updated Seurat object.
#' @importFrom Seurat CreateAssayObject DefaultAssay<-
#' @export
create_combined_assay <- function(data, combined_counts, condition_labels) {
  data[["RNA_combined"]] <- CreateAssayObject(counts = combined_counts)
  DefaultAssay(data) <- "RNA_combined"

  if (!is.null(condition_labels)) {
    data$condition <- factor(condition_labels)
  }

  timestamped_cat("Combined layers into new assay 'RNA_combined'.\n")
  return(data)
}
