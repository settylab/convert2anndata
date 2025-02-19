#' Process Metadata and Pairwise Matrices
#'
#' This function processes the metadata and pairwise matrices in a `SingleCellExperiment` object.
#' It organizes metadata into distinct categories (`uns`, `obsp`, `obsm`, `varp`, `varm`) based
#' on their dimensions and relationship to the main data matrix (`X`). Metadata matrices matching
#' observation or variable dimensions are moved into the corresponding slots, while remaining
#' metadata is retained in `uns`. Pairwise matrices for observations and variables are identified
#' and organized into `obsp` and `varp`, respectively.
#'
#' @param sce A `SingleCellExperiment` object containing the data to be processed.
#' @param alt_exps_data A list of alternative experiment data to be added to the metadata.
#' @param X The main data matrix (`anndata.X`) for the AnnData object. This is used to validate
#'   dimensions of metadata and pairwise matrices.
#' @param obsm A list of observation matrices, which can be extended during processing.
#' @return A list containing:
#'   \item{uns}{A list of remaining metadata that was not categorized into other slots.}
#'   \item{obsp}{A list of pairwise observation matrices.}
#'   \item{varp}{A list of pairwise variable matrices.}
#'   \item{obsm}{A list of observation matrices.}
#'   \item{varm}{A list of variable matrices.}
#' @importFrom SingleCellExperiment int_metadata
#' @importFrom S4Vectors metadata
#' @export
process_metadata_and_pairwise <- function(sce, alt_exps_data, X, obsm = list()) {
  timestamped_cat("Gathering metadata and pairwise matrices...\n")
  
  # Merge metadata and internal metadata
  uns <- c(metadata(sce), int_metadata(sce))
  obsp <- list()
  varp <- list()
  varm <- list()
  n_obs <- nrow(X)
  n_var <- ncol(X)
  
  # Process each metadata item
  for (name in names(uns)) {
    item <- uns[[name]]
    if (is.matrix(item) && all(dim(item) == c(n_obs, n_obs))) {
      obsp[[name]] <- item
      uns[[name]] <- NULL # Remove from uns if added to obsp
      timestamped_cat(sprintf(
        "Transferred '%s' to pairwise observations (obsp).\n",
        name
      ))
    } else if (is.matrix(item) && nrow(item) == n_obs) {
      obsm[[name]] <- item
      uns[[name]] <- NULL # Remove from uns if added to obsm
      timestamped_cat(sprintf(
        "Transferred '%s' to observation matrices (obsm).\n",
        name
      ))
    } else if (is.matrix(item) && all(dim(item) == c(n_var, n_var))) {
      varp[[name]] <- item
      uns[[name]] <- NULL # Remove from uns if added to varp
      timestamped_cat(sprintf(
        "Transferred '%s' to pairwise variables (varp).\n",
        name
      ))
    } else if (is.matrix(item) && nrow(item) == n_var) {
      varm[[name]] <- item
      uns[[name]] <- NULL # Remove from uns if added to varm
      timestamped_cat(sprintf(
        "Transferred '%s' to variable matrices (varm).\n",
        name
      ))
    } else if (is.list(item) && all(sapply(item, inherits, "SeuratCommand"))) {
      uns[[name]] <- convert_commands(item)  # Convert SeuratCommand list
      timestamped_cat(sprintf(
        "Converted SeuratCommand list '%s' for serialization.\n",
        name
      ))
    }
  }
  
  # Add alternative experiments to uns
  if (length(alt_exps_data) > 0) {
    uns[["altExperiments"]] <- alt_exps_data
  }
  
  timestamped_cat("Metadata and pairwise data organized.\n")
  
  return(list(uns = uns, obsp = obsp, varp = varp, obsm = obsm, varm = varm))
}
