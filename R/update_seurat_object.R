#' Update Seurat Object if Necessary
#'
#' Checks if the Seurat object is from an old version (v2) and updates it to the latest format.
#'
#' @param data A Seurat object.
#' @return An updated Seurat object.
#' @importFrom Seurat UpdateSeuratObject
#' @importFrom utils packageVersion
#' @export
update_seurat_object <- function(data) {
  if (!is.null(data@version) && data@version < package_version("3.0.0")) {
    timestamped_cat("Old Seurat v2 object detected. Updating to Seurat v3/v5 format...\n")
    data <- tryCatch(
      {
        Seurat::UpdateSeuratObject(data)
      },
      error = function(e) {
        stop("Failed to update Seurat v2 object: ", e$message)
      }
    )
    timestamped_cat("Seurat object successfully updated to version: ", as.character(data@version), "\n")
  }
  return(data)
}
