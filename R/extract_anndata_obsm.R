#' Extract obsm embeddings from an AnnData object
#'
#' Returns a named list of numeric matrices, one per `adata$obsm` key.
#' Each matrix has rows in `obs_names` order. Non-numeric or
#' dimension-mismatched entries are skipped with a warning.
#'
#' @param adata An AnnData object.
#' @param obs_names Character vector of cell names; used both for sanity-checking
#'   row count and for setting rownames.
#' @return A named list of numeric matrices (cells x dim).
#' @export
extract_anndata_obsm <- function(adata, obs_names = NULL) {
  obsm_keys <- anndata_mapping_keys(adata$obsm)
  if (length(obsm_keys) == 0) {
    return(list())
  }

  out <- list()
  for (k in obsm_keys) {
    mat <- tryCatch(as.matrix(adata$obsm[[k]]), error = function(e) NULL)
    if (is.null(mat)) {
      warning(sprintf("Could not convert obsm['%s'] to matrix; skipping.", k))
      next
    }
    if (!is.numeric(mat)) {
      warning(sprintf("obsm['%s'] is not numeric (type: %s); skipping.", k, typeof(mat)))
      next
    }
    if (!is.null(obs_names) && nrow(mat) != length(obs_names)) {
      warning(sprintf(
        "obsm['%s'] has %d rows but expected %d cells; skipping.",
        k, nrow(mat), length(obs_names)
      ))
      next
    }
    if (!is.null(obs_names)) rownames(mat) <- obs_names
    out[[k]] <- mat
  }
  out
}
