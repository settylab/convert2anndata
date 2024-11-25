#' Preserve Dimensional Reductions from Seurat Object
#'
#' Updates and transfers dimensional reductions from the Seurat object to the SingleCellExperiment object.
#'
#' @param data The Seurat object.
#' @param sce The SingleCellExperiment object.
#' @return The SingleCellExperiment object with updated dimensional reductions.
#' @importFrom SingleCellExperiment reducedDims<-
#' @export
preserve_reductions_seurat <- function(data, sce) {
  if (!is.null(data@reductions) && length(data@reductions) > 0) {
    timestamped_cat("Attempting to preserve dimensional reductions.\n")

    # Check if cell names match
    reduction_names <- names(data@reductions)
    for (reduction_name in reduction_names) {
      reduction <- data@reductions[[reduction_name]]
      # Update cell embeddings to match combined cells
      cells_in_reduction <- rownames(reduction@cell.embeddings)
      common_cells <- intersect(cells_in_reduction, colnames(data))
      if (length(common_cells) > 0) {
        reduction@cell.embeddings <- reduction@cell.embeddings[common_cells, , drop = FALSE]
        # Reorder to match the combined data
        reduction@cell.embeddings <- reduction@cell.embeddings[colnames(data), , drop = FALSE]
        data@reductions[[reduction_name]] <- reduction
        timestamped_cat("Updated reduction '", reduction_name, "' to match combined cells.\n")
        # Transfer to SCE
        reducedDims(sce)[[reduction_name]] <- reduction@cell.embeddings
      } else {
        timestamped_cat("No common cells found for reduction '", reduction_name, "'. Removing this reduction.\n")
        data@reductions[[reduction_name]] <- NULL # Remove the reduction if it doesn't match any cells
      }
    }
  }
  return(sce)
}
