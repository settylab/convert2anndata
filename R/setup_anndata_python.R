#' Configure the Python interpreter used by reticulate for AnnData I/O
#'
#' Call this once at the start of an R session to point reticulate at a
#' Python interpreter that has the `anndata` module installed. The
#' resolution order is:
#'
#' 1. `conda_env` (if supplied) — first tried as a conda environment name,
#'    then as a path to an environment prefix, then as a path under
#'    `~/micromamba/envs/<name>` or `~/miniconda3/envs/<name>`.
#' 2. `RETICULATE_PYTHON` environment variable — used directly as a Python
#'    binary path if it exists.
#' 3. `CONDA_PREFIX` environment variable — treated as the active
#'    conda/micromamba environment prefix.
#' 4. No-op (reticulate falls back to its own discovery).
#'
#' Failures emit a warning and continue rather than throwing, so users
#' who already have Python wired up don't have to special-case this call.
#'
#' @param conda_env Optional. A conda environment name (resolved against
#'   common locations: micromamba, miniconda3, conda) or an absolute path
#'   to an environment prefix.
#' @param required Logical; passed to the underlying reticulate calls.
#'   Defaults to `FALSE` so that misconfiguration emits a warning rather
#'   than aborting.
#' @param validate Logical. After setting the interpreter preference, run
#'   `check_anndata_python(action = "warn")` to confirm the env actually
#'   imports `anndata`. Defaults to `TRUE`. With `required = FALSE`
#'   reticulate's own validation is a no-op, so this is the only
#'   protection against silent typos in `conda_env`.
#' @return Invisibly returns a single-element character vector describing
#'   which path was used: one of `"conda_env"`, `"conda_env_resolved"`,
#'   `"reticulate_python"`, `"conda_prefix"`, or `"none"`.
#' @examples
#' \dontrun{
#' setup_anndata_python("my-anndata-env")
#' adata <- anndata::read_h5ad("pbmc.h5ad")
#' }
#' @export
setup_anndata_python <- function(conda_env = NULL, required = FALSE,
                                 validate = TRUE) {
  source_label <- "none"

  if (!is.null(conda_env) && nzchar(conda_env)) {
    source_label <- try_use_conda_env(conda_env, required = required)
  }

  if (source_label == "none") {
    python_path <- Sys.getenv("RETICULATE_PYTHON", unset = NA)
    if (!is.na(python_path) && nzchar(python_path) && file.exists(python_path)) {
      ok <- FALSE
      tryCatch({
        reticulate::use_python(python_path, required = required)
        ok <- TRUE
      }, error = function(e) {
        warning(sprintf(
          "Could not use Python at '%s': %s",
          python_path, conditionMessage(e)
        ))
      })
      if (ok) source_label <- "reticulate_python"
    }
  }

  if (source_label == "none") {
    conda_prefix <- Sys.getenv("CONDA_PREFIX", unset = NA)
    if (!is.na(conda_prefix) && nzchar(conda_prefix)) {
      ok <- FALSE
      tryCatch({
        reticulate::use_condaenv(conda_prefix, required = required)
        ok <- TRUE
      }, error = function(e) {
        warning(sprintf(
          "Could not activate CONDA_PREFIX '%s': %s",
          conda_prefix, conditionMessage(e)
        ))
      })
      if (ok) source_label <- "conda_prefix"
    }
  }

  if (isTRUE(validate) && source_label != "none") {
    suppressWarnings(check_anndata_python(action = "warn"))
  }

  invisible(source_label)
}

#' @keywords internal
try_use_conda_env <- function(conda_env, required = FALSE) {
  attempt <- function(arg) {
    tryCatch({
      reticulate::use_condaenv(arg, required = required)
      TRUE
    }, error = function(e) FALSE)
  }

  # 1. Treat as a name first.
  if (attempt(conda_env)) return("conda_env")

  # 2. Treat as a literal path.
  if (file.exists(conda_env) && attempt(conda_env)) return("conda_env")

  # 3. Resolve against common env-store locations.
  candidates <- c(
    file.path(path.expand("~/micromamba/envs"), conda_env),
    file.path(path.expand("~/miniconda3/envs"), conda_env),
    file.path(path.expand("~/anaconda3/envs"), conda_env),
    file.path(Sys.getenv("MAMBA_ROOT_PREFIX", unset = ""), "envs", conda_env)
  )
  candidates <- candidates[nzchar(candidates) & dir.exists(candidates)]
  for (cand in candidates) {
    if (attempt(cand)) return("conda_env_resolved")
  }

  warning(sprintf(
    "Could not activate conda environment '%s' as a name, a literal path, or under standard env-store locations.",
    conda_env
  ))
  "none"
}
