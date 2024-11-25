#' Attach altExps to SingleCellExperiment Object
#'
#' Attaches alternative experiments to the SingleCellExperiment object.
#'
#' @param data The Seurat object.
#' @param sce The SingleCellExperiment object.
#' @param altExp_names A vector of altExp names.
#' @return The SingleCellExperiment object with attached altExps.
#' @importFrom SingleCellExperiment altExp<-
#' @importFrom Seurat GetAssayData
#' @export
attach_alt_experiments_sce <- function(data, sce, altExp_names) {
  for (alt_name in altExp_names) {
    timestamped_cat("Attaching assay '", alt_name, "' as altExp.\n")
    # Extract the assay from Seurat object
    alt_assay <- data@assays[[alt_name]]

    # Check if counts slot/layer is empty
    if (inherits(alt_assay, "Assay5")) {
      alt_counts <- GetAssayData(alt_assay, layer = "counts")
    } else if (inherits(alt_assay, "Assay")) {
      alt_counts <- GetAssayData(alt_assay, slot = "counts")
    } else {
      next # Unsupported assay type
    }

    if (nrow(alt_counts) == 0 || ncol(alt_counts) == 0) {
      timestamped_cat("Counts is empty for assay '", alt_name, "'. Skipping altExp creation.\n")
      next
    }

    # Ensure that rownames and colnames are set
    if (is.null(rownames(alt_counts))) {
      rownames(alt_counts) <- rownames(alt_assay)
    }
    if (is.null(colnames(alt_counts))) {
      colnames(alt_counts) <- colnames(alt_assay)
    }

    # Construct colData for altExp
    alt_colData <- data@meta.data[colnames(alt_counts), , drop = FALSE]

    # Construct rowData for altExp
    alt_rowData <- data.frame(feature_id = rownames(alt_counts), row.names = rownames(alt_counts))

    # Manually construct SingleCellExperiment for the altExp
    alt_sce <- SingleCellExperiment::SingleCellExperiment(
      assays = list(counts = alt_counts),
      colData = alt_colData,
      rowData = alt_rowData
    )
    # Add annotation
    metadata(alt_sce)$annotation <- "Moved to altExp due to mismatching cells or features."
    # Attach to main sce
    SingleCellExperiment::altExp(sce, alt_name) <- alt_sce
  }
  return(sce)
}
