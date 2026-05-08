.onAttach <- function(libname, pkgname) {
  # Don't fire diagnostics in non-interactive sessions (knitting, CI, batch
  # scripts). The conversion entry points already self-check and emit a
  # tailored error if anything is wrong.
  if (!interactive()) return(invisible())
  msg <- paste0(
    "convert2anndata uses reticulate to call Python anndata. ",
    "Run convert2anndata::diagnose_anndata_python() if you hit env issues."
  )
  packageStartupMessage(msg)
}
