#' Extract obs/var names from an AnnData object
#'
#' Cell and feature names on AnnData objects are Python pandas Indexes;
#' R's `rownames()`/`colnames()` do not return them reliably. This helper
#' coerces them to character vectors via Python's `list()` builtin, falling
#' back to `names()` if that path fails.
#'
#' @param adata An AnnData object (e.g. as returned by `anndata::read_h5ad`).
#' @return A list with elements `obs_names` (cells) and `var_names` (features),
#'   each a character vector or `NULL` if extraction failed.
#' @export
anndata_names <- function(adata) {
  builtins <- reticulate::import_builtins()
  obs_names <- tryCatch(
    as.character(builtins$list(adata$obs_names)),
    error = function(e) NULL
  )
  var_names <- tryCatch(
    as.character(builtins$list(adata$var_names)),
    error = function(e) NULL
  )
  list(obs_names = obs_names, var_names = var_names)
}
