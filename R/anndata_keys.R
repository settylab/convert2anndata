#' Extract dict-like keys from an AnnData mapping (layers, obsm, uns, ...)
#'
#' AnnData's mapping attributes (`layers`, `obsm`, `uns`, ...) expose a
#' `.keys()` method. Different combinations of `reticulate` and `anndata`
#' return that result in different shapes:
#'
#' - newer reticulate auto-converts the Python `KeysView` to an R character
#'   vector (so `as.character(builtins$list(x$keys()))` runs `list()` on a
#'   char-vec and splits each name into individual characters);
#' - older versions hand back a Python object that needs `builtins$list()`
#'   to materialise.
#'
#' This helper handles both cases, plus the edge case where the mapping has
#' no `.keys()` method but `names()` works.
#'
#' @param mapping The mapping attribute (e.g. `adata$layers`).
#' @return A character vector of keys, or `character(0)` on failure.
#' @keywords internal
anndata_mapping_keys <- function(mapping) {
  raw <- tryCatch(mapping$keys(), error = function(e) NULL)

  if (is.character(raw)) {
    return(raw)
  }

  if (!is.null(raw)) {
    builtins <- reticulate::import_builtins()
    py_listed <- tryCatch(builtins$list(raw), error = function(e) NULL)
    if (is.character(py_listed)) return(py_listed)
    if (!is.null(py_listed)) {
      return(tryCatch(as.character(py_listed), error = function(e) character(0)))
    }
  }

  nm <- tryCatch(names(mapping), error = function(e) NULL)
  if (is.character(nm)) return(nm)
  character(0)
}
