#' Process Metadata and Pairwise Matrices
#'
#' This function processes the metadata and pairwise matrices in the SingleCellExperiment object.
#'
#' @param sce A SingleCellExperiment object.
#' @param alt_exps_data Data of alternative experiments.
#' @param X The main data matrix (anndata.X) for the AnnData object.
#' @return A list containing metadata and pairwise matrices.
#' @importFrom SingleCellExperiment int_metadata
#' @importFrom S4Vectors metadata
#' @export
process_metadata_and_pairwise <- function(sce, alt_exps_data, X) {
  timestamped_cat("Gathering metadata and pairwise matrices...\n")
  uns <- c(metadata(sce), int_metadata(sce))
  obsp <- list()
  varp <- list()
  varm <- list()
  n_obs <- nrow(X)
  n_var <- ncol(X)
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
    }
    if (is.matrix(item) && all(dim(item) == c(n_var, n_var))) {
      varp[[name]] <- item
      uns[[name]] <- NULL # Remove from uns if added to obsp
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
    }
  }
  if (length(alt_exps_data) > 0) {
    uns[["altExperiments"]] <- alt_exps_data
  }
  timestamped_cat("Metadata and pairwise data organized.\n")

  return(list(uns = uns, obsp = obsp, varp = varp, varm = varm))
}
