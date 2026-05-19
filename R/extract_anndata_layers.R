#' Extract layers from an AnnData object as a named list of matrices
#'
#' Each layer is converted to CsparseMatrix and given dimnames from
#' `obs_names` / `var_names`. The orientation is preserved as
#' AnnData-native (cells x features); callers that need genes x cells
#' (e.g. Seurat) should transpose.
#'
#' @param adata An AnnData object.
#' @param obs_names Character vector of cell names.
#' @param var_names Character vector of feature names.
#' @param exclude Character vector of layer names to skip (e.g. one already
#'   used as the primary `X`).
#' @return A named list of CsparseMatrix objects (cells x features). Empty
#'   list if the AnnData has no layers.
#' @importFrom methods as
#' @export
extract_anndata_layers <- function(adata, obs_names = NULL, var_names = NULL, exclude = character(0)) {
  layer_keys <- anndata_mapping_keys(adata$layers)
  layer_keys <- setdiff(layer_keys, exclude)
  if (length(layer_keys) == 0) {
    return(list())
  }

  out <- list()
  for (k in layer_keys) {
    mat <- tryCatch(adata$layers[[k]], error = function(e) NULL)
    if (is.null(mat)) {
      timestamped_cat(sprintf("WARNING: could not read layer '%s'; skipping.\n", k))
      next
    }
    mat <- ensure_csparse_matrix(mat)
    if (!methods::is(mat, "dgCMatrix")) {
      mat <- methods::as(methods::as(mat, "CsparseMatrix"), "generalMatrix")
    }
    if (!is.null(obs_names) && length(obs_names) == nrow(mat)) {
      rownames(mat) <- obs_names
    }
    if (!is.null(var_names) && length(var_names) == ncol(mat)) {
      colnames(mat) <- var_names
    }
    out[[k]] <- mat
  }
  out
}
