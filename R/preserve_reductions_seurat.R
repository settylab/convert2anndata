#' Preserve Dimensional Reductions from Seurat Object
#'
#' Updates and transfers dimensional reductions from the Seurat object to the
#' SingleCellExperiment object by subsetting and reordering the cell embeddings
#' to match the cells in the SCE. For reductions that cover only a subset of cells,
#' missing rows are padded with NA.
#'
#' @param data The Seurat object.
#' @param sce The SingleCellExperiment object.
#' @return The SingleCellExperiment object with updated dimensional reductions.
#' @importFrom SingleCellExperiment reducedDims<-
#' @export
preserve_reductions_seurat <- function(data, sce) {
  if (is.null(data@reductions) || length(data@reductions) == 0) {
    return(sce)
  }
  
  timestamped_cat("Attempting to preserve dimensional reductions.\n")
  reduction_names <- names(data@reductions)
  
  for (red_name in reduction_names) {
    red_obj <- data@reductions[[red_name]]
    
    # If no embeddings available, skip.
    if (is.null(red_obj@cell.embeddings)) {
      timestamped_cat("Reduction '", red_name, "' has no cell.embeddings. Skipping.\n")
      next
    }
    
    # Determine common cells between the SCE and the reduction.
    common_cells <- intersect(colnames(sce), rownames(red_obj@cell.embeddings))
    
    if (length(common_cells) == 0) {
      timestamped_cat("No common cells found for reduction '", red_name, "'. Skipping this reduction.\n")
      next
    }
    
    # Create a new embedding matrix: one row per cell in sce, and the same number of columns as in the original.
    new_emb <- matrix(
      NA, 
      nrow = ncol(sce), 
      ncol = ncol(red_obj@cell.embeddings),
      dimnames = list(colnames(sce), colnames(red_obj@cell.embeddings))
    )
    
    # Fill in the rows for cells that have embeddings.
    new_emb[common_cells, ] <- red_obj@cell.embeddings[common_cells, , drop = FALSE]
    
    if (length(common_cells) < length(colnames(sce))) {
      timestamped_cat("Warning: Reduction '", red_name, "' only covers ", length(common_cells),
                      " out of ", length(colnames(sce)), " cells. Padding with NA for missing cells.\n")
    }
    
    timestamped_cat("Updated reduction '", red_name, "' to match combined cells.\n")
    reducedDims(sce)[[red_name]] <- new_emb
  }
  
  return(sce)
}
