#' Extract the primary count matrix from an AnnData object
#'
#' Selects a layer from `adata$layers` if `layer` is supplied and present;
#' otherwise falls back to `adata$X`. The returned matrix is in
#' AnnData orientation (cells x features) and converted to CsparseMatrix.
#' Row and column names are set from `obs_names` / `var_names` when provided.
#'
#' @param adata An AnnData object.
#' @param layer Optional name of a layer in `adata$layers` to use as the
#'   primary matrix. If `NULL` or missing, `adata$X` is used.
#' @param obs_names Character vector of cell names (rows). Optional.
#' @param var_names Character vector of feature names (columns). Optional.
#' @return A CsparseMatrix of shape `n_obs x n_vars` (cells x features).
#' @importFrom methods as
#' @export
extract_anndata_X <- function(adata, layer = NULL, obs_names = NULL, var_names = NULL) {
  layer_keys <- anndata_mapping_keys(adata$layers)

  if (!is.null(layer) && nzchar(layer) && layer %in% layer_keys) {
    timestamped_cat(sprintf("Using layer '%s' as primary matrix.\n", layer))
    mat <- adata$layers[[layer]]
  } else {
    if (!is.null(layer) && nzchar(layer)) {
      timestamped_cat(sprintf(
        "Layer '%s' not found in adata$layers; falling back to adata$X.\n",
        layer
      ))
    } else {
      timestamped_cat("Using adata$X as primary matrix.\n")
    }
    mat <- adata$X
  }

  if (is.null(mat)) {
    stop("Could not obtain a primary matrix from the AnnData object.")
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
  mat
}
