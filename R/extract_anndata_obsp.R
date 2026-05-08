#' Extract obsp (cell x cell) matrices from an AnnData object
#'
#' AnnData stores cell-pair matrices (e.g. nearest-neighbor graphs from
#' `sc.pp.neighbors`: `connectivities`, `distances`) in `adata.obsp`.
#' This helper returns them as a named list of sparse cell x cell
#' CsparseMatrix objects, ready to be attached as Seurat `Graphs` or
#' SCE `colPairs`.
#'
#' @param adata An AnnData object.
#' @param obs_names Character vector of cell names; used both for sanity-
#'   checking shape and for setting dimnames.
#' @return A named list of sparse matrices (`n_obs x n_obs`). Empty
#'   list if `adata$obsp` is empty / missing.
#' @keywords internal
extract_anndata_obsp <- function(adata, obs_names = NULL) {
  obsp <- tryCatch(adata$obsp, error = function(e) NULL)
  if (is.null(obsp)) return(list())

  keys <- anndata_mapping_keys(obsp)
  if (length(keys) == 0) return(list())

  out <- list()
  for (k in keys) {
    mat <- tryCatch(obsp[[k]], error = function(e) NULL)
    if (is.null(mat)) {
      warning(sprintf("Could not read obsp['%s']; skipping.", k))
      next
    }
    mat <- tryCatch(ensure_csparse_matrix(mat), error = function(e) NULL)
    if (is.null(mat)) {
      warning(sprintf("Could not coerce obsp['%s'] to CsparseMatrix; skipping.", k))
      next
    }
    if (nrow(mat) != ncol(mat)) {
      warning(sprintf(
        "obsp['%s'] is not square (%d x %d); skipping.",
        k, nrow(mat), ncol(mat)
      ))
      next
    }
    if (!is.null(obs_names) && nrow(mat) != length(obs_names)) {
      warning(sprintf(
        "obsp['%s'] has %d rows but expected %d cells; skipping.",
        k, nrow(mat), length(obs_names)
      ))
      next
    }
    if (!is.null(obs_names)) {
      rownames(mat) <- obs_names
      colnames(mat) <- obs_names
    }
    out[[k]] <- mat
  }
  out
}
