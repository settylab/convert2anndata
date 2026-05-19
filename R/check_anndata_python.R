#' Check that reticulate has a working Python with the `anndata` module
#'
#' Reticulate-based packages typically fail late and cryptically when the
#' Python side isn't wired up correctly: the user gets a numpy import
#' error from deep inside Python, far from the call site, with no hint
#' as to what to fix. This helper runs a layered check at the start of
#' a conversion call and stops with a single message that says exactly
#' what is wrong and how to fix it.
#'
#' Checks performed, in order, with stop-at-first-failure semantics:
#'
#' 1. **A Python interpreter is reachable.** If `reticulate::py_available()`
#'    is `FALSE` even after initialisation, no Python could be located.
#' 2. **`PYTHONPATH` does not point at a different Python's site-packages.**
#'    A common HPC failure mode is that loading an R module exports
#'    paths like `/app/.../Python-3.11/site-packages` while reticulate
#'    has selected a 3.8 / 3.12 interpreter, producing the cryptic
#'    "Error importing numpy: you should not try to import numpy from
#'    its source directory" error. We detect this *before* it fires.
#' 3. **The Python `anndata` module is importable.**
#' 4. **`numpy` imports without complaint** (catches any remaining ABI
#'    mismatches).
#' 5. **An end-to-end smoke probe**: round-trip a 2x2 AnnData object
#'    through reticulate. This is the only way to be sure the env is
#'    actually usable for conversion work.
#'
#' Each failure produces a tailored error message including the Python
#' interpreter that was selected and a one-line suggested fix.
#'
#' @param action What to do on failure: `"stop"` (default) raises an
#'   error; `"warn"` emits a warning and returns `FALSE`; `"silent"`
#'   returns `FALSE` quietly.
#' @param verbose Logical. When `TRUE`, prints a full diagnostic block
#'   on success showing Python path, version, anndata version, numpy
#'   version, `RETICULATE_PYTHON` and `CONDA_PREFIX`. Useful as a
#'   one-stop "what does my env actually look like" call.
#' @return `TRUE` (invisibly) when all checks pass, `FALSE` when any
#'   check fails and `action != "stop"`.
#' @examples
#' \dontrun{
#' check_anndata_python()                    # stops with a helpful message
#' check_anndata_python(verbose = TRUE)      # also prints env diagnostics
#' if (check_anndata_python(action = "silent")) {
#'   adata <- anndata::read_h5ad(path)
#' }
#' }
#' @export
check_anndata_python <- function(action = c("stop", "warn", "silent"),
                                 verbose = FALSE) {
  action <- match.arg(action)

  emit <- function(msg) {
    switch(action,
      stop   = stop(msg, call. = FALSE),
      warn   = { warning(msg, call. = FALSE); FALSE },
      silent = FALSE
    )
  }

  python_path <- tryCatch(
    reticulate::py_config()$python,
    error = function(e) NA_character_
  )

  # 1. Reachable Python
  if (!isTRUE(reticulate::py_available(initialize = TRUE))) {
    return(emit(paste0(
      "No Python interpreter available to reticulate.\n",
      "Try one of:\n",
      "  Sys.setenv(RETICULATE_PYTHON = '/path/to/python')\n",
      "  reticulate::use_condaenv('your-env')\n",
      "  convert2anndata::setup_anndata_python('your-env')\n",
      "  anndata::install_anndata()  # let reticulate manage a venv"
    )))
  }

  # 2. PYTHONPATH mismatch (HPC pollution detection)
  pp_problem <- detect_pythonpath_mismatch(python_path)
  if (!is.null(pp_problem)) {
    return(emit(paste0(
      "PYTHONPATH appears to point at a different Python's site-packages.\n",
      "  Selected Python: ", python_path, "\n",
      "  Conflicting PYTHONPATH entry: ", pp_problem, "\n",
      "This is a common HPC failure mode (an R module pushes paths for\n",
      "another Python version onto PYTHONPATH). It produces the cryptic\n",
      "'Error importing numpy: you should not try to import numpy from\n",
      "its source directory' error if you proceed.\n",
      "Fix: clear PYTHONPATH before launching R:\n",
      "  unset PYTHONPATH && Rscript your-script.R\n",
      "Or from inside R, before any reticulate import:\n",
      "  Sys.unsetenv('PYTHONPATH'); .rs.restartR()  # or quit + restart"
    )))
  }

  # 3. anndata module importable
  if (!isTRUE(reticulate::py_module_available("anndata"))) {
    return(emit(paste0(
      "Python module 'anndata' is not importable.\n",
      "  Selected Python: ", python_path, "\n",
      "Fixes:\n",
      "  - Install the module into the active Python:\n",
      "      reticulate::py_install('anndata') \n",
      "  - Or use the anndata R package's bundled installer:\n",
      "      anndata::install_anndata()\n",
      "  - Or point reticulate at a Python that already has anndata:\n",
      "      Sys.setenv(RETICULATE_PYTHON = '/path/to/python')"
    )))
  }

  # 4. numpy import works
  numpy_err <- tryCatch({
    reticulate::py_run_string("import numpy")
    NULL
  }, error = function(e) conditionMessage(e))

  if (!is.null(numpy_err)) {
    msg <- paste0(
      "Python `import numpy` failed.\n",
      "  Selected Python: ", python_path, "\n",
      "  Error: ", numpy_err, "\n",
      "If the error mentions 'do not try to import numpy from its source\n",
      "directory', PYTHONPATH is likely polluted by an HPC module that\n",
      "exposes site-packages for a different Python version. Clear it:\n",
      "  Sys.unsetenv('PYTHONPATH')\n",
      "and start a fresh R session before retrying."
    )
    return(emit(msg))
  }

  # 5. End-to-end smoke probe
  smoke_err <- tryCatch({
    reticulate::py_run_string(paste(
      "import anndata, numpy as np",
      "_probe_ad = anndata.AnnData(X=np.array([[1.0,2.0],[3.0,4.0]]))",
      "_probe_shape = _probe_ad.shape",
      sep = "\n"
    ))
    NULL
  }, error = function(e) conditionMessage(e))

  if (!is.null(smoke_err)) {
    return(emit(paste0(
      "anndata smoke probe failed: a Python `anndata.AnnData(...)` call\n",
      "did not complete cleanly.\n",
      "  Selected Python: ", python_path, "\n",
      "  Error: ", smoke_err, "\n",
      "Likely causes: ABI mismatch between numpy and anndata, or a partial\n",
      "install. Reinstall into a clean env:\n",
      "  anndata::install_anndata()"
    )))
  }

  if (isTRUE(verbose)) {
    print_python_env_diagnostics(python_path)
  }

  invisible(TRUE)
}

#' Print a diagnostic block for the currently-active Python env
#'
#' Convenience wrapper for `check_anndata_python(verbose = TRUE,
#' action = "silent")`. Always prints the diagnostic block; never
#' raises. Returns `TRUE` if the env is usable, `FALSE` otherwise.
#'
#' @return Logical, invisibly.
#' @export
diagnose_anndata_python <- function() {
  ok <- suppressWarnings(check_anndata_python(action = "silent"))
  python_path <- tryCatch(
    reticulate::py_config()$python,
    error = function(e) NA_character_
  )
  print_python_env_diagnostics(python_path, ok = ok)
  invisible(ok)
}

# ---- Internals ----------------------------------------------------------

#' @keywords internal
detect_pythonpath_mismatch <- function(python_path) {
  pp <- Sys.getenv("PYTHONPATH", unset = "")
  if (!nzchar(pp) || is.na(python_path) || !nzchar(python_path)) return(NULL)

  # Determine the Python major.minor of the selected interpreter from
  # `py_config()` (most reliable) or the binary path.
  cfg_ver <- tryCatch(
    as.character(reticulate::py_config()$version),
    error = function(e) NA_character_
  )
  active_mm <- if (!is.na(cfg_ver) && nzchar(cfg_ver)) {
    sub("^(\\d+\\.\\d+).*$", "\\1", cfg_ver)
  } else NA_character_

  if (is.na(active_mm) || !nzchar(active_mm)) return(NULL)

  # Look for any entry on PYTHONPATH that pins a Python major.minor
  # (e.g. `python3.10/site-packages`, `Python-3.11/site-packages`,
  # `Python/3.11/...`) where the X.Y disagrees with the active
  # interpreter.
  patterns <- c(
    "[Pp]ython(\\d+\\.\\d+)",      # python3.10  Python3.11
    "[Pp]ython[-/](\\d+\\.\\d+)"   # Python-3.11 Python/3.11
  )
  entries <- strsplit(pp, .Platform$path.sep, fixed = TRUE)[[1]]
  for (entry in entries) {
    for (pat in patterns) {
      m <- regmatches(entry, regexec(pat, entry))[[1]]
      if (length(m) >= 2 && nzchar(m[2]) && m[2] != active_mm) {
        return(entry)
      }
    }
  }
  NULL
}

#' @keywords internal
print_python_env_diagnostics <- function(python_path, ok = TRUE) {
  cat("---- convert2anndata Python environment ----\n")
  cat("status:           ", if (isTRUE(ok)) "OK" else "BROKEN", "\n")
  cat("python:           ", python_path, "\n")
  cfg_ver <- tryCatch(reticulate::py_config()$version, error = function(e) NA)
  cat("python version:   ", as.character(cfg_ver), "\n")
  cat("RETICULATE_PYTHON:", Sys.getenv("RETICULATE_PYTHON", unset = "<unset>"), "\n")
  cat("CONDA_PREFIX:     ", Sys.getenv("CONDA_PREFIX",      unset = "<unset>"), "\n")
  pp <- Sys.getenv("PYTHONPATH", unset = "<unset>")
  cat("PYTHONPATH:       ", if (nchar(pp) > 80) paste0(substr(pp, 1, 77), "...") else pp, "\n")
  ad_v <- tryCatch(
    reticulate::py_run_string(
      "import anndata, numpy; _v = (anndata.__version__, numpy.__version__)",
      local = TRUE
    )$`_v`,
    error = function(e) NULL
  )
  if (!is.null(ad_v)) {
    cat("anndata (Python): ", ad_v[[1]], "\n")
    cat("numpy:            ", ad_v[[2]], "\n")
  } else {
    cat("anndata (Python):  <not importable>\n")
    cat("numpy:             <not importable>\n")
  }
  cat("--------------------------------------------\n")
}
