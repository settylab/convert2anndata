#' Extract `adata.raw` as a CsparseMatrix
#'
#' `adata.raw` (a snapshot of an AnnData often used by scanpy to retain
#' raw counts after normalisation) is a separate AnnData-like object
#' with its own `.X`, `.var`, and `.var_names`. It may have a different
#' number of features from the main AnnData object.
#'
#' This helper returns the raw matrix in AnnData orientation (cells x
#' features), with cell rownames from `obs_names` (the active AnnData)
#' and feature column names from `adata$raw$var_names`. Returns `NULL`
#' if `adata$raw` is unset.
#'
#' @param adata An AnnData object.
#' @param obs_names Character vector of cell names. Optional.
#' @return A CsparseMatrix of shape `n_obs x n_raw_vars`, or `NULL`.
#' @keywords internal
extract_anndata_raw <- function(adata, obs_names = NULL) {
  raw <- tryCatch(adata$raw, error = function(e) NULL)
  if (is.null(raw)) return(NULL)

  raw_X <- tryCatch(raw$X, error = function(e) NULL)
  if (is.null(raw_X)) return(NULL)

  mat <- ensure_csparse_matrix(raw_X)
  if (!methods::is(mat, "dgCMatrix")) {
    mat <- methods::as(methods::as(mat, "CsparseMatrix"), "generalMatrix")
  }

  raw_var_names <- tryCatch({
    nm <- raw$var_names
    if (is.character(nm)) nm else as.character(reticulate::import_builtins()$list(nm))
  }, error = function(e) NULL)

  if (!is.null(obs_names) && length(obs_names) == nrow(mat)) {
    rownames(mat) <- obs_names
  }
  if (!is.null(raw_var_names) && length(raw_var_names) == ncol(mat)) {
    colnames(mat) <- raw_var_names
  }
  mat
}
